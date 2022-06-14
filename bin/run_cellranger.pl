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

B<--no_seurat>              :   Do NOT run the R script to create the seurat objects.

B<--no_cleanup>             :   Do NOT remove existing files from prior runs.

B<--samplesheet_dir, -s>    :   Path to dir in which samplesheet .csv files will be written.

B<--cluster_config, -c>     :   Path to a file containing cluster resource specifications.

B<--count_reference>		: 	Specify a path to a custom reference for use with cellranger count.

B<--dry_run>                :   Create the output files but do not run.

B<--filter_single_index>    :   Pass the '--filter-single-index' argument to mkfastq

=head1 DESCRIPTION

Create and run a Snakefile that controls several cellranger actions and downstream analysis.
To wit:

    1. mkfastq
    2. count
    3. mat2csv
    4. [optional] seurat

Note that Step 4 an R script that runs the Seurat package is optional and may excluded 
with B<--no_seurat>.

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

By default, the pipeline will use the GRCh38 human reference genome for the B<--transcriptome>
argument in the cellranger count step of the pipeline.  Users may provide a path to a custom
genome reference directory to be used by the cellranger count command by passing it in with
B<--count_reference>.  For more info on how to create such references, refer to:
	https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references
and
	https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_mr

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
If it already exists, it will be used.  If B<--cluster_config> is given, the file passed 
in will be used.

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

=cut

use Cwd qw( abs_path getcwd);
use File::Copy;
use File::Path;
use FindBin qw($RealBin);
use Getopt::Long qw( :config no_auto_abbrev no_ignore_case );
use List::MoreUtils qw( first_index );
use Pod::Usage;
use Text::CSV qw( csv );

use Data::Dumper;

my $default_samplesheet_dir = "sample-csvs";  # Note: using './' for relative paths is "strongly
                                              # discouraged by snakemake for being "redundant".
my $default_cluster_config  = './cluster.config';
my $default_mapping_file    = './mapping_file';
my $seurat_cpus             = 32;
my $default_count_ref_dir   = '/fdb/cellranger/refdata-cellranger-GRCh38-1.2.0'; 

my $default_working_dir = getcwd();

my %opts;
GetOptions( \%opts,
            'mapping_file|m=s',
            'working_dir|w=s',
            'no_seurat',
            'no_cleanup',
            'samplesheet_dir|s=s',
            'cluster_config|c=s',
            'dry_run',
			'count_reference=s',
            'filter_single_index',
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

# Create cluster.config if needed.
create_cluster_config() unless ( -f $opts{ cluster_config } );

# Create Snakefile
create_snakefile( $samples, $libraries, $groups );

# Cleanup?
clean_existing_runs( $libraries );

# Copy datafiles for the R script
copy_R_datafiles() unless ( $opts{ no_seurat } );

# Run snakemake via hpc
create_snake_shell();
run_snakemake();

exit(0);


sub load_libraries {
# Read in the mapping file and store the data in a hash ref for easy access.

    my ( $libraries, @samples, $groups );
    my @csv_header = ( 'Lane', 'Sample', 'Index' );

    my @mapfile_errors;

    # Might need to remove ^M chars from the mapping file:
    if ( `grep -P '\r' $opts{ mapping_file }` ) {

        my @syscmd = ('dos2unix', $opts{ mapping_file } );
        system( @syscmd ) && die "Found ^M characters in the mapping file and couldn't remove them.\n"; 

    }

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

                for my $subgroup ( split( ',', $group ) ) {
                    push @{$groups->{ $subgroup }->{ $library }}, $sample;
                }

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

            print $cfh  "    $sample: $library"."_mkfastq/outs/fastq_path\n";

        }

    }

    print $cfh "sample_libraries:\n";
    for my $library ( keys %$libraries ) {

        for my $sample ( sort { $a cmp $b } @{$libraries->{ $library }->{ samples }} ) {

            print $cfh  "    $sample: $library\n";

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
        "time"            : "3-00:00:00"
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

    my ( $samples, $libraries, $groups ) = @_;

    my ( $target_string, $count_output, $count_dir, $seurat_output_list );

    $count_dir    = '"{sample}_count/outs/"';
    $count_output = 'directory("{sample}_count/outs/filtered_feature_bc_matrix")';

    # $target_string will vary depending on the parameters used
    my @targets;
    my @seurat_output;

    if ( $opts{no_seurat} ) {

        for my $library ( sort { $a cmp $b } keys %$libraries ) {

            for my $sample ( sort { $a cmp $b } @{$libraries->{ $library }->{ samples }} ) {

                push @targets, "'$library/outs/$sample.csv'";

            }

        }

    } else {

        # Need to add final targets for the seurat files that are to be made,
        # named for the groups and singleton samples.
        for my $group ( keys %$groups ) {

            unless ( $group eq 'NA' ) {

                for my $library ( keys %{$groups->{ $group }} ) {

                    my $pseudopath = "'$opts{ working_dir }/seurat_files/group_"."$group.merged.seurat.Rdata'";
                    push @targets, $pseudopath;
                    push @seurat_output, $pseudopath;

                }

            }

        }
        for my $library ( sort { $a cmp $b } keys %$libraries ) {

            for my $sample ( sort { $a cmp $b } @{$libraries->{ $library }->{ samples }} ) {

                my $pseudopath = "'$opts{ working_dir }/seurat_files/$sample.seurat.Rdata'";
                push @targets, $pseudopath;
                push @seurat_output, $pseudopath;

            }

        }

        $target_string = join( ",\n        ", @targets );
        $seurat_output_list = join( ",\n        ", @seurat_output);

    }


    ## This will be either the 'rule all' input or the input for the seurat step:
    my @file_list;
    for my $library ( sort { $a cmp $b } keys %$libraries ) {

        for my $sample ( sort { $a cmp $b } @{$libraries->{ $library }->{ samples }} ) {

            push @file_list, "'$library/outs/$sample.csv'";

        }

    }

    my $seurat_input_list;
    $seurat_input_list = join( ",\n        ", @file_list );

    if ( $opts{ no_seurat } ) {
        $target_string = join( ",\n        ", @file_list );
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
        directory("{library}_mkfastq/outs/fastq_path")
    shell:
        """
        module load cellranger
        if [ -d {wildcards.library}_mkfastq ]; then 
            rm -rf {wildcards.library}_mkfastq
        fi
        cellranger mkfastq --id={wildcards.library}_mkfastq \\
        --run={input.rundir} --csv={input.csv} \\);

    if ( $opts { filter_single_index } ) {
        $contents .= qq( \\
        --filter-single-index \\);
    }

    $contents .= qq(
        --localcores=12 --localmem=12
        """

rule cellranger_count:
    input:
        lambda wildcards: config["count_inputs"][wildcards.sample]
    params:
        out_dir = $count_dir,
        library = lambda wildcards: config["sample_libraries"][wildcards.sample]
    output:
        $count_output
    shell:
        """
        module load cellranger
        if [ -d {params.out_dir} ]; then rm -rf {params.out_dir}; fi
        cellranger count --id={wildcards.sample}_count \\
        --fastqs={input} --sample={wildcards.sample} \\
        --transcriptome=$opts{ count_reference } \\
        --localcores=32 --localmem=120
        mkdir -p {params.library}/outs/
        cp {wildcards.sample}_count/outs/web_summary.html \\
        {params.library}/outs/{wildcards.sample}_web_summary.html
        """

rule cellranger_mat2csv:
    input:
        "{sample}_count/outs/filtered_feature_bc_matrix"
    output:
        "{library}/outs/{sample}.csv"
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
    output:
        $seurat_output_list
    shell:
        """
        module load R
        Rscript $RealBin/CreateSeuratObjectFromDenseMatrix.R --workdir $opts{ working_dir } -c $seurat_cpus --mapfile $opts{ mapping_file }
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
snakemake --cluster-config $opts{ cluster_config } --cluster "\$sbcmd" --latency-wait=180 --jobs=32

EOF

    my $shell_file = "$opts{ working_dir }/cellranger.sh";
    open( my $sfh, '>', $shell_file ) || die "Can't write to shell script $shell_file: $!\n";

    print $sfh $contents;

}


sub copy_R_datafiles {

    my @files = ( 'Biomart_hsapiens_ensembl_gene_symbols.csv',
                  'genes_to_filter.txt',
                  'enrichr_libraries.txt' );

    my $errors = '';

    for my $filename ( @files ) {

        my $oldpath = "$RealBin/$filename";
        my $newpath = "$opts{ working_dir }/$filename";

        if ( -s $oldpath ) {

            copy $oldpath, $newpath;

        } else {

            $errors .= "Can't find $oldpath, or it is empty\n";

        }

    }

    die $errors if $errors;

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

	$opts{ count_reference } = $opts{ count_reference } // $default_count_ref_dir;

    die $errors if $errors;

}



