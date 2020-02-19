#!/usr/bin/env perl
use warnings;
use strict;
$|++;

=head1 SYNOPSIS

run_cellranger.pl - run cellranger on a set of inputs

=head1 USAGE

run_cellranger.pl -m [rundir mapping]

=head1 OPTIONS

B<--mapping_file, -m>       :   Tabular file relating libraries, samples, and groups.

B<--working_dir, -w>        :   Working directory. [Default: . (cwd)]

B<--run_aggr, -a>           :   Invoke cellranger's aggr step.

B<--no_seurat>              :   Do NOT run the R script to create the seurat objects.

B<--no_cleanup>             :   Do NOT remove existing files from prior runs.

B<--samplesheet_dir, -s>    :   Path to dir containing samplesheet .csv files per library.

B<--cluster_config, -c>     :   Path to a file containing cluster resource specifications.

B<--dry_run>                :   Create the output files but do not run.

=head1 DESCRIPTION

Create and run a Snakefile that controls several cellranger actions and downstream analysis.
To wit:

    1. mkfastq
    2. count
    3. [optional] aggr
    4. mat2csv
	5. [optional] seurat

Note that Steps 3 & 5, aggr (aggregation) and Seurat (R package) are optional.
Step 3 is not run by default, it must be explicitly included via B<--run_aggr>;
Step 5 is run by default, it must be explicitly excluded with B<--no_seurat>.
Note that these two steps are (presently) mutually exclusive, as aggregating 
the counts will create new files that are dissociated from the mapping found in
the B<--mapping_file>.

This script will look for and remove any previously existing files (snakemake 
locks and cellranger __.mro files) to allow complete reruns to take place, 
unless B<--no_cleanup> is provided.

=head1 INPUT

B<--mapping_file> should contain columns separated by tabs in the order:

    1. Library ID
    2. Sample ID
    3. Raw Data Folder
    4. Lane
    5. Sample Barcode Index
    6. Group designation

where "Library ID" is the string to be used as an identifier for the sequencing 
run found in the given "Raw Data Folder". "Sample ID" is provided to the cellranger
commands as the value for --sample.  "Lane" is what should appear in the sample-sheet's
'Lane' column, and "Sample Barcode Index" also is written to the sample sheet, as the
'Index' column value.  "Group designation" determines which samples get merged into the
same Seurat object during the Seurat step.

=head1 OUTPUT

This script will produce several files prior to execution of the workflow:

=over 1

=item B<Snakefile> - contains the commands used by snakemake to invoke the
workflow.

=item B<config.yaml> - configfile in YAML format to specify various parameters
of the workflow.  This can be re-used with the Snakefile to manually re-
instantiate the workflow without re-running this script.

=item B<cellranger.sh> - shell script that will be executed to start the pipeline.

=item B<[library].csv> - each library will get a samplesheet, based on info within
the provided B<--mapping_file>.  These .csv files will be created in B<--samplesheet_dir>

=item B<cluster.config> - file containing resource specifications for each rule. If
"cluster.config" does not exist in the current working directory, it will be created.
If it already exists, it will be used.  If B<--cluster_config> is given it will use the
file passed in.

=back

In addition to these files, of course, there will be a large amount of
data created by the cellranger command during the workflow execution.  Some of
the intermediate data will be cleaned up by this tool unless the B<--no_cleanup>
option is used.

When invoked with B<--dry_run>, the above list of files is created but the pipeline is
not executed.  This allows fine-tuning of the cluster.config paramaters, among other
things.  After editing as desired, the pipeline can be launched by running the
following command (from an sinteractive session on a biowulf node, of course):

  sbatch --cpus-per-task=2 --mem=12G --time=10-00:00:00 <work_dir>/cellranger.sh

It is important not to skimp on the time allowance for this first command.  This
command places the snakemake job on a new node, and must be able to launch and wait
for the rest of the steps to execute.  It also asks for a nominal amount of resources
to allow it to run the mat2csv cellranger command locally (that is, on the same node
as the Snakemake process).  

=head1 CONTACT

    Jason Inman
    inmanjm@nih.gov

=head1 TODO:

    1. declare intermediate files so snakemake automatically slims down the final disk usage
    2. Better fine-tuning of cluster resources [especially for the seurat step]

=cut

use Cwd qw( abs_path getcwd);
use File::Path;
use FindBin qw($RealBin);
use Getopt::Long qw( :config no_auto_abbrev no_ignore_case );
use List::MoreUtils qw( first_index );
use Pod::Usage;
use Text::CSV qw( csv );

use Data::Dumper;

my $datadir_base            = '/data/NCATS_ifx/data/Chromium';
#my $default_samplesheet_dir = "$datadir_base/sample-csvs";
my $default_samplesheet_dir = "$datadir_base/inmanjm/sample-csvs";
my $default_cluster_config  = './cluster.config';
my $default_mapping_file    = './mapping_file';
my $seurat_cpus             = 32;

my $default_working_dir = getcwd;

my %opts;
GetOptions( \%opts,
            'mapping_file|m=s',
            'working_dir|w=s',
            'run_aggr|a',
			'no_seurat',
            'no_cleanup',
            'samplesheet_dir|s=s',
            'cluster_config|c=s',
            'dry_run',
            'help|h',
           ) || die "Error getting options! $!\n";
pod2usage( -verbose => 2, -exitval => 1 ) if $opts{ help };

check_params();
chdir $opts{ working_dir } || die "$opts{ working_dir } is not accessible: $!\n";

# Load in library & sample data
my ( $libraries, $samples, $groups ) = load_libraries();

# Create config.yaml for Snakefile
create_config_file( $libraries );

# Create a samplesheet .csv for each library
create_samplesheets( $libraries );

# Create aggr's .csv if needed
my $aggr_csv = create_aggr_csv( $libraries ) if ( $opts{ run_aggr } );

# Create cluster.config if needed.
create_cluster_config() unless ( -f $opts{ cluster_config } );

# Create Snakefile
create_snakefile( $samples, $aggr_csv, $libraries, $groups );

# Cleanup?
clean_existing_runs( $libraries );

# Run snakemake via hpc
create_snake_shell();
run_snakemake();

exit(0);


sub load_libraries {
# Read in the mapping file and store the data in a hash ref for easy access.

    my ( $libraries, @samples, $groups );
    my @csv_header = ( 'Lane', 'Sample', 'Index' );

    my @mapfile_errors;

    open( my $mfh, '<', $opts{ mapping_file } ) || die "Can't open mapping file: $!\n";
    while( <$mfh> ) {

        chomp;

        my ( $library, $sample, $rundir, $lane, $index, $group ) = split( "\t", $_ );

        if ( $library && $sample && $rundir && $lane && $index ) {

            if ( not $opts{ no_seurat } and not defined($group) ) {

                push @mapfile_errors, "Found this unexpected line in mapping file (missing group?):\n$_";
                next;

            }
 
            if ( -d $rundir ) {  

                # Store rundir, begin csv:
                $libraries->{ $library }->{ rundir } = $rundir;
                unless ( exists ( $libraries->{ $library }->{ csv } ) ) {
                    push @{$libraries->{ $library }->{ csv }}, \@csv_header;
                }

            } else {

                push @mapfile_errors, "Rundir '$rundir' doesn't seem to exist/isn't a directory.";

            }

            my @csv_line = ( $lane, $sample, $index );
            push @{$libraries->{ $library }->{ csv }}, \@csv_line; 
            push @{$libraries->{ $library }->{ samples }}, $sample;
            push @samples, $sample;

            unless ( $opts{ no_seurat } ) {

                $groups->{ $group }->{ $library }++;

            }

        } else {

            push @mapfile_errors, "Found this unexpected line in mapping file:\n$_";


        }

    }

    if ( @mapfile_errors ) {

        my $epitaph = join( "\n", @mapfile_errors );
        die "Found the following errors with the mapping_file ($opts{ mapping_file }):\n$epitaph\n";

    }


    return( $libraries, \@samples, $groups );

}


sub create_config_file {
# Create config.yaml to be used by the Snakefile to easily pass in the required
# values of rundir and sample names, etc. per library

    my ( $libraries ) = @_;

    my $config_file = "$opts{ working_dir }/config.yaml";
    open( my $cfh, '>', $config_file ) || die "Can't create config_file $config_file: $!\n";
    print $cfh "rundirs:\n";
    print $cfh map{ "    $_: $libraries->{ $_ }->{ rundir }\n" } sort { $a cmp $b } keys %$libraries;

    print $cfh "count_inputs:\n";
    for my $library ( keys %$libraries ) {

        for my $sample ( sort { $a cmp $b } @{$libraries->{ $library }->{ samples }} ) {

            print $cfh  "    $sample: $library"."_mkfastq/outs/fastq_path/\n";

        }

    }

}


sub create_samplesheets {
# Create .csv sampleheets per library.  (technically, per rundir)

    my ( $libraries ) = @_;

    for my $library ( keys %$libraries ) {

        csv( in       => $libraries->{ $library }->{ csv },
             out      => "$opts{ samplesheet_dir }/$library.csv",
             sep_char => ',' );

    }

}


sub create_aggr_csv {
# Create the .csv for cellranger aggr, in the form:
#     library_id,molecule.h5
#     <library>,<path_to_library's_molecule_info.h5>

    my ( $libraries ) = @_;
    my $csv_file = "$opts{ working_dir }/aggr.csv";
    my $csv_data;
    push @$csv_data, ['library_id','molecule_h5'];

    for my $library ( keys %$libraries ) {

        for my $sample ( @{$libraries->{ $library }->{ samples }} ) {

            push @$csv_data, [ $sample, "$opts{ working_dir }/$sample"."_count/outs/molecule_info.h5" ];

        }

    }
    
    csv( in         => $csv_data,
         out        => $csv_file,
         sep_char   => ',' );

    return $csv_file;

}


sub create_cluster_config {
# Create cluster.config for specifying cluster resources
# per rule.  This allows fine-tuning of resource requests.

    my $contents = qq(
{
    "__default__" :
    {
        "cpus-per-task"   : 2,
        "mem"             : "2G",
        "time"            : "24:00:00"
    },
    "cellranger_mkfastq" :
    {
        "cpus-per-task"   : 12,
        "mem"             : "12G",
    },
    "cellranger_count" :
    {
        "cpus-per-task"   : 32,
        "mem"             : "120G",
        "time"            : "24:00:00"
    },
    "cellranger_aggr" :
    {
        "cpus-per-task"   : 12,
        "mem"             : "12G",
    },
    "R_seurat" :
    {
        "cpus-per-task"   : $seurat_cpus,
        "mem"             : "120G",
    }
}
);

    open( my $cfh, '>', $opts{ cluster_config } ) || die "Can't open $opts{ cluster_config }: $!\n";

    print $cfh $contents;

}


sub create_snakefile {
# Create the Snakefile that will be executed by snakemake.

    my ( $samples, $aggr_csv, $libraries, $groups ) = @_;

    my ( $target_string, $count_output, $count_dir, $aggr_inputs );

    $count_dir    = '"{sample}_count/outs/"';
    $count_output = 'directory("{sample}_count/outs/filtered_feature_bc_matrix")';

	# $target_string will vary depending on the parameters used
    if ( $opts{ run_aggr } ) {

        $aggr_inputs = join( ",\n        ", map{ "'$_"."_count/outs/molecule_info.h5'" } @$samples );
        $aggr_inputs =~ s/"//g;

        $count_output = '"{sample}_count/outs/molecule_info.h5"';

        $target_string = "'aggr/outs/final.csv'";

    } elsif ( $opts{no_seurat} ) {

        my @targets;
        for my $library ( sort { $a cmp $b } keys %$libraries ) {

            for my $sample ( sort { $a cmp $b } @{$libraries->{ $library }->{ samples }} ) {

                push @targets, "'$library/outs/$sample.csv'";

            }

        }
        $target_string = join( ",\n        ", @targets );

    } else {

        # Need to setup final targets as the seurat files that are to be made,
        # named for the groups 
        my @targets;
        for my $group ( keys %$groups ) {

			for my $library ( keys %{$groups->{ $group }} ) {

	            push @targets, "'$opts{ working_dir }/$library/outs/seurat_files/group_"."$group.merged.seurat.Rdata'";

			}

        }

        $target_string = join( ",\n        ", @targets );

	}

	## This will be either the 'rule all' input or the input for the seurat step:
	my @file_list;
	for my $library ( sort { $a cmp $b } keys %$libraries ) {

		for my $sample ( sort { $a cmp $b } @{$libraries->{ $library }->{ samples }} ) {

			push @file_list, "'$library/outs/$sample.csv'";

		}

	}
	my $seurat_input_list;
	if ( $opts{ no_seurat } ) {
		$target_string = join( ",\n        ", @file_list );
	} elsif ( ! $opts{ run_aggr } ) {
		$seurat_input_list = join( ",\n        ", @file_list );
	}

    my $snakefile = "$opts{ working_dir }/Snakefile"; 
    open( my $sfh, '>', $snakefile ) || die "Can't create snakefile $snakefile: $!\n";

    my $contents = qq(
configfile: "config.yaml"

localrules: cellranger_mat2csv

rule all:
    input:
        $target_string 

rule cellranger_mkfastq:
    input:
        rundir = lambda wildcards: config["rundirs"][wildcards.library],
        csv = "$opts{ samplesheet_dir }/{library}.csv"
    output:
        directory("{library}_mkfastq/outs/fastq_path/")
    shell:
        """
        module load cellranger
        rmdir -p {output}
        cellranger mkfastq --id={wildcards.library}_mkfastq \\
        --run={input.rundir} --csv={input.csv} \\
        --localcores=12 --localmem=12
        """
);

    $contents .= qq(
rule cellranger_count:
    input:
        lambda wildcards: config["count_inputs"][wildcards.sample]
    params:
        out_dir = $count_dir
    output:
        $count_output
    shell:
        """
        module load cellranger
        rmdir -p {params.out_dir}
        cellranger count --id={wildcards.sample}_count \\
        --fastqs={input} --sample={wildcards.sample} \\
        --transcriptome=/fdb/cellranger/refdata-cellranger-GRCh38-1.2.0 \\
        --localcores=32 --localmem=120
        """
);

    my ( $mat2csv_input, $mat2csv_output, $aggr_output, $aggr_remove );

    if ( $opts{ run_aggr } ) {

        $aggr_output = '"aggr/outs/filtered_feature_bc_matrix"';
        $aggr_remove = '"aggr/outs"';

        $contents .= qq(
rule cellranger_aggr:
    input:
        $aggr_inputs,
        "$aggr_csv"
    output:
        directory($aggr_output)
    shell:
        """
        module load cellranger
        rmdir -p $aggr_remove
        cellranger aggr --id=aggr \\
        --csv=$aggr_csv \\
        --normalize=mapped \\
        --localcores=12 --localmem=12
        """

);

        $mat2csv_input  = $aggr_output;
        $mat2csv_output = '"aggr/outs/final.csv"';

    } else {
        $mat2csv_input = '"{sample}_count/outs/filtered_feature_bc_matrix"';
        $mat2csv_output = '"{library}/outs/{sample}.csv"';
    }

    $contents .= qq(
rule cellranger_mat2csv:
    input:
        $mat2csv_input
    output:
        $mat2csv_output
    shell:
        """
        module load cellranger
        cellranger mat2csv {input} {output}
        """
);

	unless ( $opts{ no_seurat } ) {

		$contents .= qq(
rule R_seurat:
	input:
		$seurat_input_list
	params:
		workdir="$opts{ working_dir }/{library}/outs/"
	output:
		"$opts{ working_dir }/{library}/outs/seurat_files/group_{group}.merged.seurat.Rdata"
	shell:
		"""
		module load R
		Rscript $RealBin/CreateSeuratObjectFromDenseMatrix.R --workdir {params.workdir}	 -c $seurat_cpus --mapfile $opts{ mapping_file }
		"""
);
	}	

    # Thus, it is written.
    print $sfh $contents;

    return $snakefile;

}


sub clean_existing_runs {
# Remove *.mro files and existing directories for the given libraries

    my ( $libraries ) = @_;

    for my $library ( keys %$libraries ) {

        # Remove _count dirs and .mro's for our sample ids
        for my $sample ( @{$libraries->{ $library }->{ samples }} ) {

            my $mro_file = "__$sample" . "_count.mro";

            if ( -f $mro_file ) {

                if ( $opts{ dry_run } ) {

                    print "Found $mro_file to remove (but won't because of --dry_run).\n";

                } else {

                    print "Removing $mro_file ... ";
                    unlink $mro_file;
                    print "Done\n";

                }

            }

            my $prog_dir = "$sample" . "_count";

            if ( -d $prog_dir ) {

                if ( $opts{ dry_run } ) {

                    print "Found $prog_dir to remove (but won't because of --dry_run).\n";

                } else {

                    print "Removing $prog_dir ... ";
                    rmtree $prog_dir;
                    print "Done.\n";

                }

            }

        }

        my $mro_file = "__$library" . "_mkfastq.mro";

        if ( -f $mro_file ) {

            if ( $opts{ dry_run } ) {

                print "Found $mro_file to remove (but won't because of --dry_run).\n";

            } else {

                print "Removing $mro_file ... ";
                unlink $mro_file;
                print "Done\n";

            }

        }

        my $prog_dir = "$library" . "_mkfastq";

        if ( -d $prog_dir ) {

            if ( $opts{ dry_run } ) {

                print "Found $prog_dir to remove (but won't because of --dry_run).\n";

            } else {

                print "Removing $prog_dir ... ";
                rmtree $prog_dir;
                print "Done.\n";

            }

        }

    }

    if ( $opts{ run_aggr } ) {

        my $aggr_out = './aggr';
        if ( -d $aggr_out ) {

            if ( $opts{ dry_run } ) {

                print "Found $aggr_out to remove (but won't because of --dry_run).\n";

            } else {

                print "Removing $aggr_out ... ";
                rmtree $aggr_out;
                print "Done.\n";

            }

        }

    }

}


sub create_snake_shell {

    my $contents =<<"EOF";
#!/bin/bash

module load snakemake

# This clears any existing cellranger 'locks', left behind in the event of
# a premature cancellation of a prior run:
snakemake --unlock

## This section makes an image of the DAG for the workflow.
module load graphviz
snakemake -n --dag | dot -Tsvg > ./dag.cellranger.svg

# This is where the snakemake command goes.  See hpc documentation:
# https://hpc.nih.gov/apps/snakemake.html

sbcmd="sbatch --cpus-per-task={cluster.cpus-per-task} --mem={cluster.mem} --exclusive --constraint=x2650 --time={cluster.time}"
snakemake --cluster-config $opts{ cluster_config } --cluster "\$sbcmd" --latency-wait=120 --jobs=32

EOF

    my $shell_file = "$opts{ working_dir }/cellranger.sh";
    open( my $sfh, '>', $shell_file ) || die "Can't write to shell script $shell_file: $!\n";

    print $sfh $contents;

}


sub run_snakemake {

    my @command = ( 'sbatch', '--cpus-per-task=2', '--mem=12G', '--time=10-00:00:00' );
    push @command, "$opts{ working_dir }/cellranger.sh";

    if ( $opts{ dry_run } ) {

        print "Dry run.  Skipping execution of the pipelne but would have run:\n";
        print join( ' ', @command ) . "\n";

    } else {

        print "Now launching pipeline.\n";

        system( @command ) && die "Oops... uh: $!\n";

    }

}


sub check_params {

    my $errors = '';

    if ( $opts{ mapping_file } ) {

		$opts{ mapping_file } = abs_path( $opts{ mapping_file } );

        $errors .= "--mapping_file $opts{ mapping_file } doesn't exist or is empty!\n" 
            unless ( -s $opts{ mapping_file } );

    } else {

		$opts{ mapping_file } = abs_path( $default_mapping_file );
		if ( ! -s $opts{ mapping_file } ) {
	        $errors .= "Please provide a --mapping_file\n";
		}

    }

    $opts{ working_dir } = $opts{ working_dir } // $default_working_dir;

    $opts{ samplesheet_dir } = $opts{ samplesheet_dir } // $default_samplesheet_dir;
    mkdir $opts{ samplesheet_dir } unless ( -d $opts{ samplesheet_dir } );

    $opts{ cluster_config } = $opts{ cluster_config } // $default_cluster_config;

    die $errors if $errors;

}



