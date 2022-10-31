#

# Author: Naveed Ishaque
# Date 31/01/2016
# Requirements : 2 cores, ~4Gb mem

# PRE - input is path to results_per_pid directory
# POST - produces PID and gene recurrency matrix based of default SNV and INDEL folder and file names
 
########################################
# Version 0.1
# Date: 01/02/2016
# Performs basic gene counts
# Very thin documentation
########################################
# Version 0.2
# Date: 01/10/2016
# Makes onco print table
# Very thin documentation
########################################
# Version 0.3
# Date: 05/10/2016
# USES A FUNCTION!!!!!!!
########################################
# Version 0.4
# Date: 05/10/2016
# Bug fix - previously was displaying intronic mutations :/
########################################
# Version 0.5
# Date: 06/10/2016
# Identifies SHM, ncRNA, uses MMML colours
########################################
# Version 0.6
# Date: 07/10/2016
# use getopts
########################################
# Version 0.7-K20K
# Date: 20/10/2017
# Modified for K20K and including SVs
########################################
# Version 0.10
# Date 19 Feb 2018
# Changed default kataegis to 6 per kb (as per Alexandrov 2013)
# parse arm level cnvs
# parse gene level cnvs
# modified the reporting of SVs
# Removed "simple" mode
# print to predefined output file
# print to predefined sample_info file
# added exclusion list
########################################
# Version 11
# Date Feb 2018
# Added more stuff ?
# filtered homodels with no SNPs
########################################
# Version 12
# Date April 2018
# Intogen tested and works
#######################################
# Version 13
# Date April 2018
# Parse sample level information (can now deal with multi sample pids)
# Removed multi named genes
########################################
# Version 14
# Date June 21 2018
# Modify how CNAs are parsed based on meeting on T3601
########################################
# Version 15
# Date July 2018
# Includes rare germline variants (snvs and indels)
# Date December 2018
# LOH events are now disentangled
# Date March 2019
# Genes with duplicate names are now not skipped, but printed to a separate oncoprint table
########################################

# TODO: accept tabluar input
# TODO: add SV type
# TODO: refine identification of focal amplifcation
# TODO: refine kataegis file (region + annotate genes)
# TODO: create kataegis plots
# TODO: report everything! Expect that R-script will filter
# TODO: count number of events for each mutation type

use strict;
use List::MoreUtils qw(uniq);
use Getopt::Long;
use Cwd qw(cwd);

my $usage = "\nThis script produces PID and gene recurrency matrix based of default MPILEUP-SNV, PLATYPUS-INDEL, and SOPHIA-SV folder and file names\n\n\tperl $0 -r [/path/to/results_per_pid] -m [min gene recurrence]\n\n";
my $version = "v0.15";
my $version_noDots = "v015"; # else intogen will not work

#######################
## PROJECT VARIABLES ##
#######################

my $results_dir = "/icgc/dkfzlsdf/analysis/hipo/my_hipo_project/whole_genome_sequencing/results_per_pid/";
my $genome_fai = "/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/sequence/1KGRef/hs37d5.fa.fai";
my $chrom_arms = "/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/stats/hg19armsizesXY.plain.bed";
my $genes_bed = "/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/databases/gencode/gencode19/gencode.v19.annotation_plain.genes.bed";

####################
## FILE VARIABLES ##
####################

# extensions / suffix (use either the default OTP files or any of the edited files listed below):
# snvs and indels:
my $ext_snv = "_somatic_snvs_conf_8_to_10.vcf"; # comment: default from OTP 
#my $ext_snv = "_somatic_functional_snvs_conf_8_to_10.vcf"; # comment: use only if you do not care about kataegis
#my $ext_snv = "_somatic_snvs_conf_8_to_10_blacklistRemoved_TiN_rescued.vcf"; # comment: from BODA/CO blacklist/TiNDA rescue workflow
#my $ext_snv = "_somatic-TiN_snvs_conf_8_to_10.artifactsRemoved.vcf";
#my $ext_snv = "_somatic_functional_snvs_conf_8_to_10_blacklistRemoved_TiN_rescued.vcf"; # comment: from BODA blacklist/TiNDA rescue workflow, use only if you do not care about kataegis 

my $ext_indel = "_somatic_indels_conf_8_to_10.vcf"; # comment: default from OTP
#my $ext_indel = "_somatic_functional_indels_conf_8_to_10.vcf"; # comment: use only if you do not care about kataegis
#my $ext_indel = "_somatic_indels_conf_8_to_10_blacklistRemoved_TiN_rescued.vcf"; # comment: from BODA/CO blacklist/TiNDA rescue workflow 
#my $ext_indel = "_somatic-TiN_indels_conf_8_to_10.artifactsRemoved.head.vcf";
#my $ext_indel = "_somatic_functional_indels_conf_8_to_10_blacklistRemoved_TiN_rescued.vcf"; # comment: from BODA blacklist/TiNDA rescue workflow, use only if you do not care about kataegis 


# germline (use the script TiNDA_filterHighConfiGermlineRareVariants.pl	to create the input files):
my $ext_snv_germline = "_germline_functional_TiN_filtered_conf_8_to_10.vcf";
my $ext_indel_germline = "_germline_functional_TiN_filtered_regi_conf_8_to_10.vcf";


# svs and cnvs (use either the default OTP files or any of the edited files listed below):
my $ext_sv  = "_filtered_somatic_minEventScore3.tsv"; # comment: default from OTP 
#my $ext_sv = "_*_filtered_somatic.tsv"; # comment: default from OTP, lower confidence
my $ext_cnv = "_most_important_info*"; # comment: default from OTP


# prefix
my $pre_snv   = "snvs_";
my $pre_indel = "indel_";
my $pre_sv    = "svs_";
my $pre_cnv   = "";

# folders
my $folder_snv   = "mpileup_";
my $folder_indel = "platypus_indel_";
my $folder_sv    = "SOPHIA_";
my $folder_cnv   = "ACEseq_";
my $folder_arriba = "fusions_arriba";

# parameters
my $shm_number = 6;
my $min_recurrence = 1;
#my $top_cnv_recurrence = 6; # not used anywhere else
my $sv_close="100000";
my $chr_arm_level_CNV_percentage = 0.3; 
my $threshold_cnv = 0.3;

# exclusion list
my %exclusion_list=();
$exclusion_list{"U3"}=1;
$exclusion_list{"U8"}=1;
$exclusion_list{"snoU13"}=1;
$exclusion_list{"Y_RNA"}=1;
$exclusion_list{"Metazoa_SRP"}=1;
$exclusion_list{"SNORD112"}=1;
$exclusion_list{"SNORA40"}=1;
$exclusion_list{"SNORA51"}=1;
$exclusion_list{"7SK"}=1;
$exclusion_list{"SNORD74"}=1;
$exclusion_list{"5SrRNA"}=1;
$exclusion_list{"ACA64"}=1;
$exclusion_list{"SCARNA11"}=1;
$exclusion_list{"SCARNA15"}=1;
$exclusion_list{"SNORA77"}=1;
$exclusion_list{"U1"}=1;
$exclusion_list{"5S_rRNA"}=1;

# technically not to be edited beyond this point =========================================================

###################
## PARSE OPTIONS ##
###################

my $pid_list;
my @pids;
my $output_prefix;
my $cnv_overview_table;
my $rnaseq_rpp_dir;

GetOptions("shm=i"                 => \$shm_number,
           "min_recurrence=i"      => \$min_recurrence,
           "kataegis=i"            => \$shm_number,
           "sv_close=i"            => \$sv_close,
           "arm_percentage=f"      => \$chr_arm_level_CNV_percentage,
           "results_dir=s"         => \$results_dir,
           "pids=s"                => \$pid_list,
           "output_prefix=s"       => \$output_prefix,
           "cnv_overview_table=s"  => \$cnv_overview_table,
           "fusion_rnaseq_dir=s"   => \$rnaseq_rpp_dir,
           "threshold_cnv=f"       => \$threshold_cnv,
           "genome_fai=s"            => \$genome_fai,
           "chrom_arms=s"            => \$chrom_arms,
           "genes_bed=s"             => \$genes_bed
           ) or die ("Error in command line arguments\n");

warn "WARNING: output_prefix not defined\n" unless defined $output_prefix;
warn "WARNING: could not find results directory: $results_dir\n"  unless -e $results_dir;
warn "WARNING: could not find genome fai file: $genome_fai\n"     unless -e $genome_fai;
warn "WARNING: could not find chromosome arm file: $chrom_arms\n" unless -e $chrom_arms;
warn "WARNING: could not find genes bed file: $genes_bed\n"       unless -e $genes_bed;

exit unless defined $output_prefix;
exit unless -e $results_dir;
exit unless -e $genome_fai;
exit unless -e $chrom_arms;
exit unless -e $genes_bed;

#@pids = split (" ",$pid_list) if (length $pid_list > 1);
#@pids = `ls $results_dir` unless (length $pid_list > 1);
my @pids;

if (length $pid_list > 1){
  warn "WARNING: variable \$pid_list was defined as \"$pid_list\"...\n";
  if (-e $pid_list){
       warn "WARNING: variable \$pid_list was defined as \"$pid_list\"... which is a file... using PID definition in file\n";
       open(my $pid_fh, "<$pid_list") or die "ERROR: cannot open file $pid_list";	
       while(my $pid_line = <$pid_fh>){
         chomp $pid_line;
         my @pids_line_arr = split (" ",$pid_line);
         push(@pids, @pids_line_arr);
       }
       close($pid_fh);
       warn "\#\# PIDs found in $pid_list: ".join(" ",@pids)."\n";
  }
  else {
       warn "WARNING: variable \$pid_list was defined as \"$pid_list\"... which is not a file... parsing space separated PIDs\n";
       @pids = split (" ",$pid_list);
  }
}
else{
  @pids = `ls $results_dir`
}


chomp @pids;
my $num_pids = scalar @pids;

my $outfile             = "$output_prefix.$version.min$min_recurrence.kataegis$shm_number.sv$sv_close.cnv$chr_arm_level_CNV_percentage.onco_print.tsv";
my $outfile_dupl        = "$output_prefix.$version.min$min_recurrence.kataegis$shm_number.sv$sv_close.cnv$chr_arm_level_CNV_percentage.onco_print_dupl.tsv";
my $outfile_sample_info = "$output_prefix.$version.min$min_recurrence.kataegis$shm_number.sv$sv_close.cnv$chr_arm_level_CNV_percentage.sample_info.tsv";
my $outfile_kataegis    = "$output_prefix.$version.min$min_recurrence.kataegis$shm_number.sv$sv_close.cnv$chr_arm_level_CNV_percentage.kataegis.tsv";
my $intogen_dir         = "${output_prefix}_${version_noDots}_min${min_recurrence}_kataegis${shm_number}_sv${sv_close}_kataegis_INTOGEN";
my $files_used          = "$output_prefix.$version.min$min_recurrence.kataegis$shm_number.sv$sv_close.cnv$chr_arm_level_CNV_percentage.files_used.tsv";
my $log_file            = "$output_prefix.$version.min$min_recurrence.kataegis$shm_number.sv$sv_close.cnv$chr_arm_level_CNV_percentage.log.txt";

##########################
## print the parameters ##
##########################

warn "\n\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n\n";
warn "\# PROJECT DIR AND PIDS \#\n";
warn "GENOME_FAI:     $genome_fai\n";
warn "CHROM_ARMS:     $chrom_arms\n";
warn "GENES_BED:      $genes_bed\n";
warn "RESULTS_DIR:    $results_dir\n";
warn "PIDS:           @pids\n";
warn "NUM_PIDS:       $num_pids\n";
warn "\n";
warn "\# FOLDER/FILE WILDCARDS \#\n";
warn "SNV FILE:       PID/$folder_snv/$pre_snv"."[PID]"."$ext_snv\n";
warn "SNV FILE GERMLINE:       PID/$folder_snv/$pre_snv"."[PID]"."$ext_snv_germline\n";
warn "INDEL FILE:     PID/$folder_indel/$pre_indel"."[PID]"."$ext_indel\n";
warn "INDEL FILE GERMLINE:     PID/$folder_indel/$pre_indel"."[PID]"."$ext_indel_germline\n";
warn "SV FILE:        PID/$folder_sv/$pre_sv"."[PID]"."$ext_sv\n";
warn "CNV FILE:       PID/$folder_cnv/$pre_cnv"."[PID]"."$ext_cnv\n";
warn "\n";
warn "\# VARIABLES \#\n";
warn "MIN_RECURRENCE: $min_recurrence\n";
warn "KATAEGIS:       $shm_number (per 1kb)\n";
warn "SV_CLOSE:       $sv_close\n";
warn "ARM_LEVEL_CNV   ".($chr_arm_level_CNV_percentage*100)."% of arm\n";
warn "\n";
warn "\# OUTFILES \#\n";
warn "GENERECURRENCE: $outfile\n";
warn "DUPLGENERECURRENCE: $outfile_dupl\n";
warn "SAMPLE INFO:    $outfile_sample_info\n";
warn "KATAEGIS:       $outfile_kataegis\n";
warn "INTOGEN DIR:    $intogen_dir\n";
warn "RNASEQ_DIR:     $rnaseq_rpp_dir\n";
warn "FILES USED:     $files_used\n";
warn "LOG FILE:       $log_file\n";
warn "\n\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n\n";

###############################
## initiate output variables ##
###############################

my %output_hash;
my %output_hash_germline;
#my %output_hash_pid; # not used anywhere
my %sample_info;
my %kataegis;

#########################
## GLOBAL RESULTS VARS ##
#########################
my %mutation_files;

#########################
## Check files per PID ##
#########################

warn "\# CHECKING: PIDS, SAMPLES + SNV, INDEL, SV and CNV files ...\n";

my @pid_samples_temp;

open (my $files_used_fh, ">$files_used") or die "\n\nERROR: cannot open file: $files_used\n\n";
open (my $log_file_fh, ">$log_file")     or die "\n\nERROR: cannot open file: $log_file\n\n";


## Reading cnv overview table
my %cnv_overview = ();
if($cnv_overview_table !~/^$/) {
  open(CNV_OVERVIEW, "$cnv_overview_table") || die "can't open the $cnv_overview_table\n";
  while(<CNV_OVERVIEW>) {
    chomp;
    if($_!~/^PID/) {
      my @ss = split(/\t/, $_);
      $cnv_overview{$ss[0]}{'ploidy_factor'} = $ss[2];
      $cnv_overview{$ss[0]}{'purity'} = $ss[3];
    }
  }
  close CNV_OVERVIEW;
}

### Looping over the pids
foreach my $pid (@pids){
  chomp $pid;

  my @mpileup_dirs = `ls $results_dir/$pid/$folder_snv* 2>/dev/null`;
  if (length($mpileup_dirs[0]) < 3 ){
    warn "WARNING: cannot find mpileup dir for $pid ($results_dir/$pid/$folder_snv\*) .... skipping\n";
    next;
  }
  
  chomp(my @sample_combinations = `ls -d $results_dir/$pid/$folder_snv* | xargs -n1 basename`);

  my @sample_combinations_new = map {s/$folder_snv//r; } @sample_combinations;

  foreach my $sample_combination(@sample_combinations_new) {

    print "Processing: $pid\t$sample_combination";
    print $log_file_fh "\nProcessing: $pid\t$sample_combination\n";

    my $pid_sample = $pid."_".$sample_combination;
    push (@pid_samples_temp, $pid_sample);

    my %expected_files;
    $expected_files{"mutation_file"}="$results_dir/$pid/$folder_snv$sample_combination/$pre_snv$pid$ext_snv";
    $expected_files{"mutation_file_germline"}="$results_dir/$pid/$folder_snv$sample_combination/$pre_snv$pid$ext_snv_germline";
    $expected_files{"indel_file"}="$results_dir/$pid/$folder_indel$sample_combination/$pre_indel$pid$ext_indel";
    $expected_files{"indel_file_germline"}="$results_dir/$pid/$folder_indel$sample_combination/$pre_indel$pid$ext_indel_germline";
    $expected_files{"sv_file"}="$results_dir/$pid/$folder_sv$sample_combination/$pre_sv$pid$ext_sv";
    $expected_files{"cnv_file"}="$results_dir/$pid/$folder_cnv$sample_combination/$pre_cnv$pid$ext_cnv";

    $mutation_files{$pid_sample}{"pid"}                    = "$pid";
    $mutation_files{$pid_sample}{"sample"}                 = "$sample_combination";
    $mutation_files{$pid_sample}{"mutation_file"}          = (`ls $expected_files{"mutation_file"} 2>/dev/null`);#          if (-e $expected_files{"mutation_file"});
    $mutation_files{$pid_sample}{"snv_file"}               = (`ls $expected_files{"mutation_file"} 2>/dev/null`);#          if (-e $expected_files{"mutation_file"});
    $mutation_files{$pid_sample}{"mutation_file_germline"} = (`ls $expected_files{"mutation_file_germline"} 2>/dev/null`);# if (-e $expected_files{"mutation_file_germline"});
    $mutation_files{$pid_sample}{"snv_file_germline"}      = (`ls $expected_files{"mutation_file_germline"} 2>/dev/null`);# if (-e $expected_files{"mutation_file_germline"});
    $mutation_files{$pid_sample}{"indel_file"}             = (`ls $expected_files{"indel_file"} 2>/dev/null`);#             if (-e $expected_files{"indel_file"});
    $mutation_files{$pid_sample}{"indel_file_germline"}    = (`ls $expected_files{"indel_file_germline"} 2>/dev/null`);#    if (-e $expected_files{"indel_file_germline"});
    $mutation_files{$pid_sample}{"sv_file"}                = (`ls $expected_files{"sv_file"} 2>/dev/null`);#                if (-e $expected_files{"sv_file"});

### CNV option
    if($cnv_overview_table=~/^$/) {
      my @ACEseq_most_importants                             = (`ls $expected_files{"cnv_file"} 2>/dev/null`);#               if ( scalar( () = <dir/$expected_files{"cnv_file"}> ) );
      if (scalar (@ACEseq_most_importants) > 1) {
#       warn "WARNING: found more than one ACEseq solution for $pid_sample... using lowest ploidy solution\n";
        print $log_file_fh "-- WARNING: found more than one ACEseq solution for $pid_sample... using lowest ploidy solution\n";
      }
      $mutation_files{$pid_sample}{"cnv_file"}               = $ACEseq_most_importants[0];
    }
    else {      
      $mutation_files{$pid_sample}{"cnv_file"}               = (`ls $results_dir/$pid/$folder_cnv*$sample_combination*/$pre_cnv$pid$ext_cnv$cnv_overview{$pid_sample}{'ploidy_factor'}_$cnv_overview{$pid_sample}{'purity'}.txt 2>/dev/null`); 
    }

    ### RNAseq fusions
    my @sample_combination_split = ();
    if($rnaseq_rpp_dir !~/^$/) {
      # extracting first part of the sample_combination
      @sample_combination_split = split(/_/, $sample_combination);
      $expected_files{"rnaseq_fusion_file"} = "$rnaseq_rpp_dir/$pid/$folder_arriba/$sample_combination_split[0]_${pid}.fusions.txt";
      $mutation_files{$pid_sample}{"rnaseq_fusion_file"} = (`ls $expected_files{"rnaseq_fusion_file"} 2>/dev/null`);
      
    }
    chomp $mutation_files{$pid_sample}{"mutation_file"};
    chomp $mutation_files{$pid_sample}{"mutation_file_germline"};
    chomp $mutation_files{$pid_sample}{"indel_file"};
    chomp $mutation_files{$pid_sample}{"indel_file_germline"};
    chomp $mutation_files{$pid_sample}{"sv_file"};
    chomp $mutation_files{$pid_sample}{"cnv_file"};
    chomp $mutation_files{$pid_sample}{"rnaseq_fusion_file"};

    print " [no SNV]"            unless -e $mutation_files{$pid_sample}{"mutation_file"};
    print " [no SNV germline]"   unless -e $mutation_files{$pid_sample}{"mutation_file_germline"};
    print " [no indel]"          unless -e $mutation_files{$pid_sample}{"indel_file"};
    print " [no indel germline]" unless -e $mutation_files{$pid_sample}{"indel_file_germline"};
    print " [no SV]"             unless -e $mutation_files{$pid_sample}{"sv_file"};
    print " [no CNV]"            unless -e $mutation_files{$pid_sample}{"cnv_file"};
    print "\n";

    print $log_file_fh "-- WARNING: could not find SNV file: ".           $expected_files{"mutation_file"}."\n"          unless -e $mutation_files{$pid_sample}{"mutation_file"};
    print $log_file_fh "-- WARNING: could not find SNV germline file: ".  $expected_files{"mutation_file_germline"}."\n" unless -e $mutation_files{$pid_sample}{"mutation_file_germline"};
    print $log_file_fh "-- WARNING: could not find INDEL file: ".         $expected_files{"indel_file"}."\n"             unless -e $mutation_files{$pid_sample}{"indel_file"};
    print $log_file_fh "-- WARNING: could not find INDEL germline file: ".$expected_files{"indel_file_germline"}."\n"    unless -e $mutation_files{$pid_sample}{"indel_file_germline"};
    print $log_file_fh "-- WARNING: could not find SV file: ".            $expected_files{"sv_file"}."\n"                unless -e $mutation_files{$pid_sample}{"sv_file"};
    print $log_file_fh "-- WARNING: could not find CNV file: ".           $expected_files{"cnv_file"}."\n"               unless -e $mutation_files{$pid_sample}{"cnv_file"};

    print $files_used_fh "$pid\t$sample_combination\t$pid_sample\tSNV\t".$mutation_files{$pid_sample}{"mutation_file"}."\n";
    print $files_used_fh "$pid\t$sample_combination\t$pid_sample\tSNV_germline\t".$mutation_files{$pid_sample}{"mutation_file_germline"}."\n";
    print $files_used_fh "$pid\t$sample_combination\t$pid_sample\tINDEL\t".$mutation_files{$pid_sample}{"indel_file"}."\n";
    print $files_used_fh "$pid\t$sample_combination\t$pid_sample\tINDEL_germline\t".$mutation_files{$pid_sample}{"indel_file_germline"}."\n";
    print $files_used_fh "$pid\t$sample_combination\t$pid_sample\tSV\t".$mutation_files{$pid_sample}{"sv_file"}."\n";
    print $files_used_fh "$pid\t$sample_combination\t$pid_sample\tCNV\t".$mutation_files{$pid_sample}{"cnv_file"}."\n";
    print $files_used_fh "$pid\t$sample_combination_split[0]\t$pid\tRNAseq_fusion\t".$mutation_files{$pid_sample}{"rnaseq_fusion_file"}."\n";

  }

}

close ($files_used_fh);

warn "\n\#\# Found ".scalar(@pid_samples_temp)." samples from ".scalar(@pids)." pids... switching pid-sample identifiers...\n\n";
print $log_file_fh "\n\n-- Found ".scalar(@pid_samples_temp)." samples from ".scalar(@pids)." pids... switching pid-sample identifiers...\n";

@pids = @pid_samples_temp;

warn "\# DONE! Found all files\n\n";
print $log_file_fh  "-- DONE! Found all files\n\n";

close($log_file_fh);

##################
## PARSE GENOME ##
##################

warn "\# PARSING: chromosome arms file: $chrom_arms\n";
my %chr_arm_sizes = ();

open (my $chr_arms_fh, "cat $chrom_arms |" ) or die "# ERROR: cannot open chromosome arm file: $chrom_arms\n";
  while(<$chr_arms_fh>){
    chomp;
    if (/^(.*?)\t(.*?)\t(.*?)\t(.*?)$/){
      $chr_arm_sizes{$4}=$3-$2;
    }
  }
close ($chr_arms_fh);
warn "\# DONE!\n\n";

#################
## PARSE GENES ##
#################

warn "\# PARSING: gene file: $genes_bed\n";
my %genes = ();

open (my $genes_fh, "intersectBed -wa -wb -a $genes_bed -b $chrom_arms |") or die "# ERROR: cannot open genes BED file: $genes_bed\n";
  while(<$genes_fh>){
    chomp;
    if (/^(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)$/){
      $genes{$4}{"ENSG"}  =$genes{$4}{"ENSG"}.$7;
      $genes{$4}{"chr"}   =$1 ;
      $genes{$4}{"start"} =$2 ;
      $genes{$4}{"end"}   =$3 ;
      $genes{$4}{"chrarm"}=$12;
      $genes{$4}{"count"} =$genes{$4}{"count"} + 1;
    }
  }
close ($genes_fh);

warn "\# DONE!\n";

warn "\n\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n\n";


###################
## PARSE INTOGEN ##
###################
warn "\# PARSING: intogen command (run separately) ...\n";

`mkdir -p $intogen_dir/`;
`mkdir -p $intogen_dir/vcfs`;
`mkdir -p $intogen_dir/config`;
`mkdir -p $intogen_dir/config/tmp`;
`mkdir -p $intogen_dir/MutSig`;
`mkdir -p $intogen_dir/tmp`;

my $cwd = cwd();

my $scheduler_conf = "# Default maximum number of simultaneous jobs
jobs = 8

# Default cluster queues
#queue = short-high, normal, long

[oncodrivefm_genes]

# Number of cores to use when running this task
cores = 8

# You can run some tasks in different queues
#queue = short-high, normal, long

[oncodrivefm_pathways]

# Number of cores to use when running this task
cores = 8

# You can run some tasks in different queues
#queue = short-high, normal, long
";

my $system_conf = 'temporal_dir = "'.$cwd.'/'.$intogen_dir.'/tmp"
perl_path = "perl"
perl_lib = "%(bgdata://intogen/software/1.0)/vep"
pool_analysis_enabled = True
pool_name = "all_projects"
vep_path = "%(bgdata://intogen/software/1.0)/vep/variant_effect_predictor.pl"
vep_cache = "%(bgdata://datasets/vepcache/70)"
ma_db = "%(bgdata://intogen/ma/1.0)"
transfic_path = "%(bgdata://intogen/transfic/1.0)"
kegg_path = "%(bgdata://intogen/kegg/1.0)/ensg_kegg.tsv"
gene_transcripts_path = "%(bgdata://intogen/ensembl/1.0)/gene_transcripts.tsv"
symbol_mapping_file = "%(bgdata://intogen/ensembl/1.0)/genes.tsv"
known_drivers_file = "%(bgdata://intogen/ensembl/1.0)/known_drivers.tsv"
network_file = "%(bgdata://intogen/networks/1.0)"
consequence_rankings = "%(bgdata://intogen/ensembl/1.0)/consequence_ranking.tsv"
liftover_path = "%(bgdata://datasets/liftover/hg18_to_hg19)"
mutsig_enabled = True
mutsig_path = "'.$cwd.'/'.$intogen_dir.'/MutSig"
mutsig_sym2gene_file = "%(bgdata://intogen/ensembl/1.0)/hgnc_mutsigsym_mapping.tsv"
mutsig_vep_dictionary = "%(bgdata://intogen/ensembl/1.0)/consequence_classification_dictionary.txt"
mutsig_exome_coverage = "%(bgdata://intogen/mutsigcov/1.0)/exome_full192.coverage.txt"
mutsig_covariates = "%(bgdata://intogen/mutsigcov/1.0)/gene.covariates.txt"
mutsig_hg19_path = "%(bgdata://datasets/genomereference/hg19)"
matlab_mcr = "/ibios/tbi_cluster/13.1/x86_64/matlab/MATLAB_Compiler_Runtime_R2013a/v81"
';

my $task_conf ='[annotation]
assembly = "hg19"

[oncodrivefm]
significance_threshold = 0.1
samples_threshold = 2
genes_filter_file = "${bgdata_filters}/gene-filter.default.txt.cgc"

[oncodriveclust]
significance_threshold = 0.1
samples_threshold = 5
genes_filter_file = "${bgdata_filters}/gene-filter.default.txt.cgc"

[mutsig] 
significance_threshold = 0.1
genes_filter_file = "${bgdata_filters}/gene-filter.default.txt.cgc"
';

my $task_validation="[annotation]
assembly = option('hg19')

[oncodrivefm] 
significance_threshold = float
samples_threshold = integer
genes_filter_file = file(default=None)

[oncodriveclust]
significance_threshold = float
samples_threshold = integer
genes_filter_file = file(default=None)

[mutsig]
significance_threshold = float
genes_filter_file = file(default=None)

[__many__]

  [[oncodrivefm]]
  significance_threshold = float(default=None)
  samples_threshold = integer(default=None)
  genes_filter_file = file(default=None)

  [[oncodriveclust]]
  significance_threshold = float(default=None)
  samples_threshold = integer(default=None)
  genes_filter_file = file(default=None)

  [[mutsig]]
  significance_threshold = float(default=None)
  genes_filter_file = file(default=None)

";

open (my $task_validation_fh, ">$intogen_dir/config/task.validation" ) or die "\# ERROR: cannot open file to write: $intogen_dir/conf/task.validation";
open (my $task_conf_fh,       ">$intogen_dir/config/task.conf" ) or die "\# ERROR: cannot open file to write: $intogen_dir/conf/task.conf";
open (my $scheduler_conf_fh,  ">$intogen_dir/config/scheduler.conf" ) or die "\# ERROR: cannot open file to write: $intogen_dir/conf/scheduler.conf";
open (my $system_conf_fh,     ">$intogen_dir/config/system.conf" ) or die "\# ERROR: cannot open file to write: $intogen_dir/conf/system.conf";

print $task_validation_fh "$task_validation\n";
print $task_conf_fh       "$task_conf\n";
print $scheduler_conf_fh  "$scheduler_conf\n";
print $system_conf_fh     "$system_conf\n";

close($task_validation_fh);
close($task_conf_fh);
close($scheduler_conf_fh);
close($system_conf_fh);

`cp -R /icgc/ngs_share/general/intogen/MutSigCV_1.4/* $intogen_dir/MutSig && chmod 777 -R $intogen_dir/MutSig`;

my $vcf_file_list = "";

foreach my $pid (@pids){
  my $snv_file =   $mutation_files{$pid}{"mutation_file"};
  my $indel_file = $mutation_files{$pid}{"indel_file"};
  chomp $snv_file;
  chomp $indel_file;

  `grep -v "^\#" $snv_file  $indel_file | sed 's|:|\t|' | cut -f 2-6 | sort -k1g -k2n | grep -v ^hs37d5 | grep -v ^GL > $intogen_dir/vcfs/$pid.snv_indel.vcf`;
  $vcf_file_list = $vcf_file_list."," if (length $vcf_file_list > 5);
  $vcf_file_list = $vcf_file_list."vcfs/$pid.snv_indel.vcf";
}

my $smconfig="id=MyProject
files=$vcf_file_list

[annotation]
assembly = \"hg19\"

[oncodrivefm]
significance_threshold = 0.1
samples_threshold = 2
genes_filter_file = \"\${bgdata_filters}/gene-filter.default.txt.cgc\"

[oncodriveclust]
significance_threshold = 0.1
samples_threshold = 5
genes_filter_file = \"\${bgdata_filters}/gene-filter.default.txt.cgc\"

[mutsig]
significance_threshold = 0.1
genes_filter_file = \"\${bgdata_filters}/gene-filter.default.txt.cgc\"";

open (my $smconfig_fh,     ">$intogen_dir/MyProject.smconfig" ) or die "\# ERROR: cannot open file to write: $intogen_dir/MyProject.smconfig";
print $smconfig_fh "$smconfig\n";
close($smconfig_fh);

my $intogen_command = "qsub -I -X -l nodes=1:ppn=8 -l mem=20g

cd $cwd/$intogen_dir

module load python/3.4.3

export INTOGEN_HOME=$cwd/$intogen_dir/config

export DSERVICES_FSLOCK_VERBOSE=1
export DSERVICES_FSLOCK_WARNING=1
export MCR_INHIBIT_CTF_LOCK=1

intogen -i MyProject.smconfig --verbose --split-size 5000";

open (my $command_fh,     ">$intogen_dir/MyProject.batch" ) or die "\# ERROR: cannot open file to write: $intogen_dir/MyProject.batch";
print $command_fh "$intogen_command\n";
close($command_fh);

warn "\# DONE!\n";
warn "\n\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n\n";

################
## PARSE CNVS ##
################

warn "\# PARSING: arm and gene level CNAs ...\n";

# TODO: identify arm level CNV/LOH per sample

my %arm_level_cnvs= ();
my %gene_level_cnvs= ();

foreach my $pid (@pids){
  my @header;
  my $chr_arm_intersect_header_chromosome = "\#chromosome";
  my $chr_arm_intersect_header_start = "start";
  my $chr_arm_intersect_header_end = "end";
  my $chr_arm_intersect_header_TCN = "TCN";
  my $chr_arm_intersect_header_length = "length";
  my $chr_arm_intersect_header_covRatio = "covRatio";
  my $chr_arm_intersect_header_C1 = "c1Mean";
  my $chr_arm_intersect_header_C2 = "c2Mean";
  my $chr_arm_intersect_header_genotype = "genotype";
  #my $chr_arm_intersect_header_GNL = "GNL";
  my $chr_arm_intersect_header_CNAtype = "CNA.type";
  my $chr_arm_intersect_header_NbrOfHetsSNPs = "NbrOfHetsSNPs";
  my %header_index;

  my $cnv_file = $mutation_files{$pid}{"cnv_file"};

  ########### TCC PLOIDY SEX #############

  my $sex = `grep "assumed sex" $cnv_file | cut -f 2 -d ":"`;
  chomp $sex;
  my $roundPloidy = `grep "roundPloidy" $cnv_file | cut -f 2 -d ":"`;
  chomp $roundPloidy;
  my $tcc = `grep "tcc" $cnv_file | cut -f 2 -d ":"`;
  chomp $tcc;
  
  $sample_info{"CNA TCC"}{$pid}=$tcc;
  $sample_info{"CNA ploidy"}{$pid}=$roundPloidy;
  $sample_info{"CNA sex"}{$pid}=$sex;

  ########### ARM #############

  open (my $chr_arm_cnv_intersect, "intersectBed -header -wa -wb -a $cnv_file -b $chrom_arms | ") or die "\# ERROR: failed intersectBed\n";
  while (<$chr_arm_cnv_intersect>){
    if (/^\#chromosome/){
      # parse header line
      chomp;
      @header = split ("\t", $_);
      my $header_counter = 0;
      foreach (@header){
        $header_index{$_}=$header_counter;
        $header_counter++;
      }
    }
    if (/^\#/){
      next;
    }
    else{
      chomp;
      my @line = split ("\t",$_);
      my $this_chr_arm = $line[(scalar @line)-1];
      my $s_cnatype = $line[$header_index{$chr_arm_intersect_header_CNAtype}];
      my $s_cnaTCN = $line[$header_index{$chr_arm_intersect_header_TCN}];
      my $s_ploidy = $sample_info{"CNA ploidy"}{$pid};
      my $s_C1 = $line[$header_index{$chr_arm_intersect_header_C1}];
      my $s_numHetSNPS = $line[$header_index{$chr_arm_intersect_header_NbrOfHetsSNPs}];

      # use thresholds suggested in T3601, and then adapted in T3869
      my $isHomoDel = ($s_cnaTCN < 0.5 && !($s_cnaTCN eq "NA")); # && !($s_C1 eq "NA"));
      my $isLOH = ($s_C1 < 0.3 && !($s_C1 eq "NA"));
      my $isLoss = (($s_ploidy - $s_cnaTCN) > $threshold_cnv && !($s_ploidy eq "NA") && !($s_cnaTCN eq "NA")); # && !($s_C1 eq "NA"));
      my $isGain = (($s_cnaTCN - $s_ploidy) > $threshold_cnv && !($s_cnaTCN eq "NA") && !($s_ploidy eq "NA")); # && !($s_C1 eq "NA"));

      # if ChrX or Y, use ACEseq CNAtype column
      if($this_chr_arm =~/^(chr)?X[pq]|^(chr)?Y[pq]/i){
        $isHomoDel = ($s_cnatype =~ "HomoDel");
        $isLOH = ($s_cnatype =~ "LOH");
        $isLoss = ($s_cnatype =~ "DEL");
        $isGain = ($s_cnatype =~ "DUP");
      }

      if ($isGain){
        # print "$pid\tgain\t$line[$header_index{$chr_arm_intersect_header_length}]\t$line[$header_index{$chr_arm_intersect_header_covRatio}]\t$this_chr_arm\n";
        $arm_level_cnvs{$this_chr_arm}{"GAIN"}{$pid} = $arm_level_cnvs{$this_chr_arm}{"GAIN"}{$pid} + $line[$header_index{$chr_arm_intersect_header_length}];
      }
      if($isHomoDel && $s_numHetSNPS >= 2){
        $arm_level_cnvs{$this_chr_arm}{"HOMODEL"}{$pid} = $arm_level_cnvs{$this_chr_arm}{"HOMODEL"}{$pid} + $line[$header_index{$chr_arm_intersect_header_length}];
      }
      if($isLoss && !($isHomoDel)){
        # print "$pid\tloss\t$line[$header_index{$chr_arm_intersect_header_length}]\t$line[$header_index{$chr_arm_intersect_header_covRatio}]\t$this_chr_arm\n";
        $arm_level_cnvs{$this_chr_arm}{"LOSS"}{$pid} = $arm_level_cnvs{$this_chr_arm}{"LOSS"}{$pid} + $line[$header_index{$chr_arm_intersect_header_length}];
      }
      if($isLOH && !($isHomoDel)){
        # print "$pid\tLOH\t$line[$header_index{$chr_arm_intersect_header_length}]\t$line[$header_index{$chr_arm_intersect_header_GNL}]\t$this_chr_arm\n";
        $arm_level_cnvs{$this_chr_arm}{"LOH"}{$pid} = $arm_level_cnvs{$this_chr_arm}{"LOH"}{$pid} + $line[$header_index{$chr_arm_intersect_header_length}];
      }
    }
  }
  close ($chr_arm_cnv_intersect);

  ############ GENES ##############

  open (my $gene_cnv_intersect, "intersectBed -header -wa -wb -a $cnv_file -b $genes_bed |") or die "\# ERROR: failed intersectBed\n";
  while (<$gene_cnv_intersect>){
    if (/^\#chromosome/){
      # parse header line
      chomp;
      @header = split ("\t", $_);
      my $header_counter = 0;
      foreach (@header){
        $header_index{$_}=$header_counter;
        $header_counter++;
      }
    }
    if (/^\#/){
      next;
    }
    else{
      chomp;
      my @line = split ("\t",$_);
      my $this_gene = $line[(scalar @line)-5]; # this is so dangerous
      my $s_cnatype = $line[$header_index{$chr_arm_intersect_header_CNAtype}];
      my $s_cnaTCN = $line[$header_index{$chr_arm_intersect_header_TCN}];
      my $s_ploidy = $sample_info{"CNA ploidy"}{$pid};
      my $s_C1 = $line[$header_index{$chr_arm_intersect_header_C1}];
      my $s_numHetSNPS = $line[$header_index{$chr_arm_intersect_header_NbrOfHetsSNPs}];

      my $isHomoDel = ($s_cnaTCN < 0.5 && !($s_cnaTCN eq "NA")); # && !($s_C1 eq "NA"));
      my $isLOH = ($s_C1 < 0.3 && !($s_C1 eq "NA"));
      my $isLoss = (($s_ploidy - $s_cnaTCN) > $threshold_cnv && !($s_ploidy eq "NA") && !($s_cnaTCN eq "NA")); # && !($s_C1 eq "NA"));
      my $isGain = (($s_cnaTCN - $s_ploidy) > $threshold_cnv && !($s_cnaTCN eq "NA") && !($s_ploidy eq "NA")); # && !($s_C1 eq "NA"));
      my $isHighGain = ($s_cnaTCN > ($s_ploidy * 2.5) && !($s_cnaTCN eq "NA") && !($s_ploidy eq "NA")); # && !($s_C1 eq "NA"));

      # if ChrX or Y, use ACEseq CNAtype column
      if($line[0] =~/^(chr)?X|^(chr)?Y/i){
        $isHomoDel = ($s_cnatype =~ "HomoDel");
        $isLOH = ($s_cnatype =~ "LOH");
        $isLoss = ($s_cnatype =~ "DEL");
        $isGain = ($s_cnatype =~ "DUP");
      }

      if($isHighGain) {
        $gene_level_cnvs{$this_gene}{"HIGHGAIN"}{$pid} = $gene_level_cnvs{$this_gene}{"HIGHGAIN"}{$pid} + 1;
      }
      if($isGain && !($isHighGain)) {
        $gene_level_cnvs{$this_gene}{"GAIN"}{$pid} = $gene_level_cnvs{$this_gene}{"GAIN"}{$pid} + 1;
      }
      if($isHomoDel && $s_numHetSNPS >= 2) {
        $gene_level_cnvs{$this_gene}{"HOMODEL"}{$pid} = $gene_level_cnvs{$this_gene}{"HOMODEL"}{$pid} + 1;
      }
      if($isLoss && !($isHomoDel)) {
        $gene_level_cnvs{$this_gene}{"LOSS"}{$pid} = $gene_level_cnvs{$this_gene}{"LOSS"}{$pid} + 1;
      }
      if($isLOH && !($isHomoDel)) {
        # LOH
        $gene_level_cnvs{$this_gene}{"LOH"}{$pid} = $gene_level_cnvs{$this_gene}{"LOH"}{$pid} + 1;
      }
    }
  }
  close ($gene_cnv_intersect);
  
  ##### TODO: identify cytoband level CNV/LOH per sample

}

#########################
## PARSE CNVS - review ##
#########################

foreach my $this_arm (sort {$a <=> $b} keys %arm_level_cnvs){
  my %cnv_list;
  $cnv_list{"GAIN"} = 0;
  $cnv_list{"LOSS"} = 0;
  $cnv_list{"LOH"} = 0;
  $cnv_list{"HOMODEL"} = 0;
  $cnv_list{"HIGHGAIN"} = 0;

  foreach my $this_pid (sort @pids){
    $cnv_list{"GAIN"}= $cnv_list{"GAIN"} + 1 if (($arm_level_cnvs{$this_arm}{"GAIN"}{$this_pid})/($chr_arm_sizes{$this_arm}) >=  $chr_arm_level_CNV_percentage);
    $cnv_list{"HIGHGAIN"}= $cnv_list{"HIGHGAIN"} + 1 if (($arm_level_cnvs{$this_arm}{"HIGHGAIN"}{$this_pid})/($chr_arm_sizes{$this_arm}) >=  $chr_arm_level_CNV_percentage);
    $cnv_list{"LOSS"}= $cnv_list{"LOSS"} + 1 if (($arm_level_cnvs{$this_arm}{"LOSS"}{$this_pid})/($chr_arm_sizes{$this_arm}) >=  $chr_arm_level_CNV_percentage);
    $cnv_list{"HOMODEL"}= $cnv_list{"HOMODEL"}+1 if (($arm_level_cnvs{$this_arm}{"HOMODEL"}{$this_pid})/($chr_arm_sizes{$this_arm}) >=  $chr_arm_level_CNV_percentage);
    $cnv_list{"LOH"} = $cnv_list{"LOH"}  + 1 if (($arm_level_cnvs{$this_arm}{"LOH"}{$this_pid}) /($chr_arm_sizes{$this_arm}) >=  $chr_arm_level_CNV_percentage);


    $sample_info{"CNA chr$this_arm"}{$this_pid} = $sample_info{"CNA chr$this_arm"}{$this_pid}."amp;" if (($arm_level_cnvs{$this_arm}{"GAIN"}{$this_pid})/($chr_arm_sizes{$this_arm}) >=  $chr_arm_level_CNV_percentage);
    $sample_info{"CNA chr$this_arm"}{$this_pid} = $sample_info{"CNA chr$this_arm"}{$this_pid}."highAmp;" if (($arm_level_cnvs{$this_arm}{"HIGHGAIN"}{$this_pid})/($chr_arm_sizes{$this_arm}) >=  $chr_arm_level_CNV_percentage);
    $sample_info{"CNA chr$this_arm"}{$this_pid} = $sample_info{"CNA chr$this_arm"}{$this_pid}."del;"  if (($arm_level_cnvs{$this_arm}{"LOSS"}{$this_pid})/($chr_arm_sizes{$this_arm}) >=  $chr_arm_level_CNV_percentage);
    $sample_info{"CNA chr$this_arm"}{$this_pid} = $sample_info{"CNA chr$this_arm"}{$this_pid}."homoDel;"  if (($arm_level_cnvs{$this_arm}{"HOMODEL"}{$this_pid})/($chr_arm_sizes{$this_arm}) >=  $chr_arm_level_CNV_percentage);
    $sample_info{"CNA chr$this_arm"}{$this_pid} = $sample_info{"CNA chr$this_arm"}{$this_pid}."LOH;"  if (($arm_level_cnvs{$this_arm}{"LOH"}{$this_pid}) /($chr_arm_sizes{$this_arm}) >=  $chr_arm_level_CNV_percentage);
  }
  warn "\# CNVs - $this_arm\tgain:".$cnv_list{"GAIN"}."\thighGain:".$cnv_list{"HIGHGAIN"}."\tloss:".$cnv_list{"LOSS"}."\thomoDel:".$cnv_list{"HOMODEL"}."\tLOH:".$cnv_list{"LOH"}."\n";
}

foreach my $this_gene (sort {$a <=> $b} keys %gene_level_cnvs){
  my %cnv_list;
  $cnv_list{"GAIN"} = 0;
  $cnv_list{"LOSS"} = 0;
  $cnv_list{"LOH"} = 0;
  $cnv_list{"HOMODEL"} = 0;
  $cnv_list{"HIGHGAIN"} = 0;

  my $this_arm = $genes{$this_gene}{"chrarm"};

  foreach my $this_pid (@pids){
  #  $cnv_list{"GAIN"}= $cnv_list{"GAIN"} + 1 if ( ($gene_level_cnvs{$this_gene}{"GAIN"}{$this_pid} eq 1) && ($arm_level_cnvs{$this_arm}{"GAIN"}{$this_pid})/($chr_arm_sizes{$this_arm}) <  $chr_arm_level_CNV_percentage);
  #  $cnv_list{"LOSS"}= $cnv_list{"LOSS"} + 1 if ( ($gene_level_cnvs{$this_gene}{"LOSS"}{$this_pid} eq 1) && ($arm_level_cnvs{$this_arm}{"LOSS"}{$this_pid})/($chr_arm_sizes{$this_arm}) <  $chr_arm_level_CNV_percentage);
  #  $cnv_list{"LOH"} = $cnv_list{"LOH"}  + 1 if ( ($gene_level_cnvs{$this_gene}{"LOH"}{$this_pid} eq 1)  && ($arm_level_cnvs{$this_arm}{"LOH"}{$this_pid}) /($chr_arm_sizes{$this_arm}) <  $chr_arm_level_CNV_percentage);

    setup_gene_snv($this_gene, $this_pid) if (!(defined $output_hash{$this_gene}));

    my $isHomoDel = $gene_level_cnvs{$this_gene}{"HOMODEL"}{$this_pid} ge 1;
    my $isLOH = $gene_level_cnvs{$this_gene}{"LOH"}{$this_pid} ge 1;
    my $isLoss = $gene_level_cnvs{$this_gene}{"LOSS"}{$this_pid} ge 1;
    my $isGain = $gene_level_cnvs{$this_gene}{"GAIN"}{$this_pid} ge 1;
    my $isHighGain = $gene_level_cnvs{$this_gene}{"HIGHGAIN"}{$this_pid} ge 1;

    my $isChrGain = (($arm_level_cnvs{$this_arm}{"GAIN"}{$this_pid})/($chr_arm_sizes{$this_arm}) >= $chr_arm_level_CNV_percentage) ;
    my $isChrLoss = (($arm_level_cnvs{$this_arm}{"LOSS"}{$this_pid})/($chr_arm_sizes{$this_arm}) >= $chr_arm_level_CNV_percentage) ;
    my $isChrLOH = (($arm_level_cnvs{$this_arm}{"LOH"}{$this_pid})/($chr_arm_sizes{$this_arm}) >= $chr_arm_level_CNV_percentage) ;
    my $isChrHomoDel = (($arm_level_cnvs{$this_arm}{"HOMODEL"}{$this_pid})/($chr_arm_sizes{$this_arm}) >= $chr_arm_level_CNV_percentage) ;

    $output_hash{$this_gene}{"chrGAIN"}{"."}{$this_pid} = $output_hash{$this_gene}{"chrGAIN"}{"."}{$this_pid} +1  if ( $isGain && $isChrGain );
    $output_hash{$this_gene}{"chrLOSS"}{"."}{$this_pid} = $output_hash{$this_gene}{"chrLOSS"}{"."}{$this_pid} +1  if ( $isLoss  && $isChrLoss);
    $output_hash{$this_gene}{"chrHOMODEL"}{"."}{$this_pid} = $output_hash{$this_gene}{"chrHOMODEL"}{"."}{$this_pid} +1  if ( $isHomoDel  && $isChrHomoDel);
    $output_hash{$this_gene}{"chrLOH"}{"."}{$this_pid}  = $output_hash{$this_gene}{"chrLOH"}{"."}{$this_pid}  +1  if ( $isLOH && $isChrLOH);

    # See Phabricator T3880	
    $output_hash{$this_gene}{"HIGHGAIN"}{"."}{$this_pid} = $output_hash{$this_gene}{"HIGHGAIN"}{"."}{$this_pid} +1  if ($isHighGain);
    #$output_hash{$this_gene}{"GAIN"}{"."}{$this_pid} = $output_hash{$this_gene}{"GAIN"}{"."}{$this_pid} +1  if ( $isGain && !$isChrGain );
    $output_hash{$this_gene}{"GAIN"}{"."}{$this_pid} = $output_hash{$this_gene}{"GAIN"}{"."}{$this_pid} +1  if ( $isGain ); # && !$isChrGain );
    #$output_hash{$this_gene}{"LOSS"}{"."}{$this_pid} = $output_hash{$this_gene}{"LOSS"}{"."}{$this_pid} +1  if ( $isLoss && !$isChrLoss);
    $output_hash{$this_gene}{"LOSS"}{"."}{$this_pid} = $output_hash{$this_gene}{"LOSS"}{"."}{$this_pid} +1  if ( $isLoss ); # && !$isChrLoss);
    #$output_hash{$this_gene}{"HOMODEL"}{"."}{$this_pid} = $output_hash{$this_gene}{"HOMODEL"}{"."}{$this_pid} +1  if ( $isHomoDel && !$isChrHomoDel);
    $output_hash{$this_gene}{"HOMODEL"}{"."}{$this_pid} = $output_hash{$this_gene}{"HOMODEL"}{"."}{$this_pid} +1  if ( $isHomoDel ); # && !$isChrHomoDel);
    #$output_hash{$this_gene}{"LOH"}{"."}{$this_pid}  = $output_hash{$this_gene}{"LOH"}{"."}{$this_pid}  +1  if ( $isLOH && !$isChrLOH);
    $output_hash{$this_gene}{"LOH"}{"."}{$this_pid}  = $output_hash{$this_gene}{"LOH"}{"."}{$this_pid}  +1  if ( $isLOH ); # && !$isChrLOH);

  }
  #warn "\# - $this_gene\tgain:".$cnv_list{"GAIN"}."\tloss:".$cnv_list{"LOSS"}."\tLOH:".$cnv_list{"LOH"}."\n";
}

warn "\# DONE!\n";
warn "\n\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n\n";



###################
## PARSE SNV NUM ##
###################

warn "\# PARSING: number of SNVs ...\n";

my $snvs_ANNOVAR_FUNCTION = "ANNOVAR_FUNCTION";
my $snvs_GENE = "GENE";
my $snvs_EXONIC_CLASSIFICATION = "EXONIC_CLASSIFICATION";

foreach my $pid (@pids){
  my $mutation_file = $mutation_files{$pid}{"mutation_file"};
  my $mutation_file_germline = $mutation_files{$pid}{"mutation_file_germline"};

  my @snvs_loci;
  my %snvs_header_index;

  my $snvs_header_line=`head -1000 $mutation_file | grep "^\#CHROM"`;
  my @snvs_header = split ("\t", $snvs_header_line);
  my $snvs_header_elem_counter=0;
  foreach (@snvs_header){
    $snvs_header_index{$_}=$snvs_header_elem_counter;
    $snvs_header_elem_counter++;
  }

  my @snvs_loci_germline;
  my %snvs_header_index_germline;

  if (-e $mutation_files{$pid}{"mutation_file_germline"}){  
    # technically, the headers should be the same for somatic and germline files, but just to be sure:
    my $snvs_header_line_germline=`head -1000 $mutation_file_germline | grep "^\#CHROM"`;
    my @snvs_header_germline = split ("\t", $snvs_header_line_germline);
    my $snvs_header_elem_counter_germline=0;
    foreach (@snvs_header_germline){
      $snvs_header_index_germline{$_}=$snvs_header_elem_counter_germline;
      $snvs_header_elem_counter_germline++;
    }
  }

  #@snvs_loci = `grep -v "^\#" $mutation_file `;

  my $snv_summary_command = "grep -v $snvs_ANNOVAR_FUNCTION $mutation_file | cut -f ".($snvs_header_index{$snvs_ANNOVAR_FUNCTION}+1).",".($snvs_header_index{$snvs_EXONIC_CLASSIFICATION}+1). " | sort | uniq -c | sed 's|^ *||g' | sed 's| |\t|'";
  my @snv_summary = `$snv_summary_command`;
  
  my @snv_summary_germline;

  if (-e $mutation_files{$pid}{"mutation_file_germline"}){
    my $snv_summary_command_germline = "grep -v $snvs_ANNOVAR_FUNCTION $mutation_file_germline | cut -f ".($snvs_header_index_germline{$snvs_ANNOVAR_FUNCTION}+1).",".($snvs_header_index_germline{$snvs_EXONIC_CLASSIFICATION}+1). " | sort | uniq -c | sed 's|^ *||g' | sed 's| |\t|'";
    @snv_summary_germline = `$snv_summary_command_germline`;
  }

  my $total_count = 0;
  foreach (@snv_summary){
    chomp;
    my @snv_locus = split ("\t", $_);
    $snv_locus[1] = "intergenic" if ($snv_locus[1] =~ m/stream/);
    $snv_locus[1] = $snv_locus[1]." - ".$snv_locus[2]  if ($snv_locus[1] =~ m/^exonic/);
    $sample_info{"SNV - $snv_locus[1]"}{$pid}=$sample_info{"SNV - $snv_locus[1]"}{$pid}+$snv_locus[0];
    $total_count = $total_count + $snv_locus[0];
  }

  my $total_count_germline = 0;
  foreach (@snv_summary_germline){
    chomp;
    my @snv_locus_germline = split ("\t", $_);
    $snv_locus_germline[1] = "intergenic" if ($snv_locus_germline[1] =~ m/stream/);
    $snv_locus_germline[1] = $snv_locus_germline[1]." - ".$snv_locus_germline[2]  if ($snv_locus_germline[1] =~ m/^exonic/);
    $sample_info{"SNV_germline - $snv_locus_germline[1]"}{$pid}=$sample_info{"SNV_germline - $snv_locus_germline[1]"}{$pid}+$snv_locus_germline[0];
    $total_count_germline = $total_count_germline + $snv_locus_germline[0];
  }
  $sample_info{"SNV - total"}{$pid}=$total_count;
  $sample_info{"SNV_germline - total"}{$pid}=$total_count_germline;
}

warn "\# DONE!\n";
warn "\n\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n\n";

#####################
## PARSE INDEL NUM ##
#####################

warn "\# PARSING: number of indels ...\n";

my $indels_ANNOVAR_FUNCTION = "ANNOVAR_FUNCTION";
my $indels_GENE = "GENE";
my $indels_EXONIC_CLASSIFICATION = "EXONIC_CLASSIFICATION";

foreach my $pid (@pids){

  my $indel_file = $mutation_files{$pid}{"indel_file"};
  my $indel_file_germline = $mutation_files{$pid}{"indel_file_germline"};

  my @indels_loci;
  my %indels_header_index=();

  my $indels_header_line=`head -1000 $indel_file | grep "^\#CHROM"`;
  my @indels_header = split ("\t", $indels_header_line);
  my $indels_header_elem_counter=0;
  foreach (@indels_header){
    $indels_header_index{$_}=$indels_header_elem_counter;
    $indels_header_elem_counter++;
  }	

  # technically, the headers should be the same for somatic and germline files, but just to be sure:
  my @indels_loci_germline;
  my %indels_header_index_germline=();

  if (-e $mutation_files{$pid}{"indel_file_germline"}){
    my $indels_header_line_germline=`head -1000 $indel_file_germline | grep "^\#CHROM"`;
    my @indels_header_germline = split ("\t", $indels_header_line_germline);
    my $indels_header_elem_counter_germline=0;
    foreach (@indels_header_germline){
      $indels_header_index_germline{$_}=$indels_header_elem_counter_germline;
      $indels_header_elem_counter_germline++;
    }
  }

  #@indels_loci = `grep -v "^\#" $indel_file `;

  my $indel_summary_command = "grep -v $indels_ANNOVAR_FUNCTION $indel_file | cut -f ".($indels_header_index{$indels_ANNOVAR_FUNCTION}+1).",".($indels_header_index{$indels_EXONIC_CLASSIFICATION}+1). " | sort | uniq -c | sed 's|^ *||g' | sed 's| |\t|'";
  my @indel_summary = `$indel_summary_command`;

  my @indel_summary_germline;

  if (-e $mutation_files{$pid}{"indel_file_germline"}){
    my $indel_summary_command_germline = "grep -v $indels_ANNOVAR_FUNCTION $indel_file_germline | cut -f ".($indels_header_index_germline{$indels_ANNOVAR_FUNCTION}+1).",".($indels_header_index_germline{$indels_EXONIC_CLASSIFICATION}+1). " | sort | uniq -c | sed 's|^ *||g' | sed 's| |\t|'";
    @indel_summary_germline = `$indel_summary_command_germline`;
  }
  
  my $total_count = 0;
  foreach (@indel_summary){
    chomp;
    my @indel_locus = split ("\t", $_);
    $indel_locus[1] = "intergenic" if ($indel_locus[1] =~ m/stream/);
    $indel_locus[1] = $indel_locus[1]." - ".$indel_locus[2]  if ($indel_locus[1] =~ m/^exonic/);
    $sample_info{"INDEL - $indel_locus[1]"}{$pid}=$sample_info{"INDEL - $indel_locus[1]"}{$pid}+$indel_locus[0];
    $total_count = $total_count + $indel_locus[0];
  }
  my $total_count_germline = 0;
  foreach (@indel_summary_germline){
    chomp;
    my @indel_locus_germline = split ("\t", $_);
    $indel_locus_germline[1] = "intergenic" if ($indel_locus_germline[1] =~ m/stream/);
    $indel_locus_germline[1] = $indel_locus_germline[1]." - ".$indel_locus_germline[2]  if ($indel_locus_germline[1] =~ m/^exonic/);
    $sample_info{"INDEL_germline - $indel_locus_germline[1]"}{$pid}=$sample_info{"INDEL_germline - $indel_locus_germline[1]"}{$pid}+$indel_locus_germline[0];
    $total_count_germline = $total_count_germline + $indel_locus_germline[0];
  }
  $sample_info{"INDEL - total"}{$pid}=$total_count;
  $sample_info{"INDEL_germline - total"}{$pid}=$total_count_germline;
}

warn "\# DONE!\n";
warn "\n\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n\n";

########################
## PARSE KATAEGIS/SHM ##
########################

warn "\# PARSING: sample specific kataegis loci ...\n";

my %shm;

foreach my $pid (@pids){

  $shm{$pid};

  my $mutation_file = $mutation_files{$pid}{"mutation_file"};
  my $indel_file = $mutation_files{$pid}{"indel_file"};

  my @header;
  my %header_index=();
  my $gene_header="GENE";
  my @shm_loci;

  @shm_loci = `grep -vw "\#" $mutation_file | cut -f 1-2 | awk -F'\\t' '{print \$1"\\t"(\$2-1)"\\t"(\$2+1000)}' | bedtools sort -i - |  bedtools intersect -wa -header -c -b - -a $mutation_file`; 

  foreach (@shm_loci){
    chomp;
    if (/^#CHROM/){
      @header = split("\t", $_);
      my $header_elem_counter=0;
      foreach my $elem (@header){
        $header_index{$elem}=$header_elem_counter;
        $header_elem_counter++;
      } 
    }
    else{
      my @line = split ("\t", $_);
      next if ($line[(scalar @line)-1] < $shm_number);
      my $locus = $line[$header_index{$gene_header}];
      if (($locus =~ tr/\;//) > 0){   
        my @snv1 = split ("\;",$locus);
        for (my $i=0; $i <= length(@snv1); $i++){
          my @genes = split (",", @snv1[$i]);
          foreach my $gene (@genes){
            next if $gene =~ m/^ENST00/; # skip if reporting multiple transcripts per gene
            if (($gene =~ tr/\(//) > 0){
              $gene =~ m/^(.*?)\(.*$/; 
              $gene = $1;
            }
            # inititate_hash
            setup_gene_snv($gene, $pid) if (!(defined $output_hash{$gene}));
            $output_hash{$gene}{"SHM"}{"."}{$pid} = "1";
#            warn "\# FOUND SHM: $pid : $gene\n";
          }
        }
      }
      # if no ";" then check for ",";
      else {
      my @genes = split (",", $locus);
      foreach my $gene (@genes){
        next if $gene =~ m/^ENST00/; # skip if reporting multiple transcripts per gene
          if (($gene =~ tr/\(//) > 0){ 
            $gene =~ m/^(.*?)\(.*$/;  
            $gene = $1;
          }
      
          # inititate_hash
          setup_gene_snv($gene, $pid) if (!(defined $output_hash{$gene}));
          $output_hash{$gene}{"SHM"}{"."}{$pid} = "1";
#          warn "\# FOUND SHM: $pid : $gene\n";
        }
      } 
    } 
  }
}
warn "\# DONE!\n";
warn "\n\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n\n";


###############
## PARSE SVS ##
###############

warn "\# PARSING: SOPHIA SVs\n";

my %svs;
my $svs_gene1="gene1"; # normally column 19 or index 18
my $svs_gene2="gene2"; # normally column 29 or index 28
my $svs_nearestCodingGeneUpstream1     = "nearestCodingGeneUpstream1";
my $svs_nearestCancerGeneUpstream1     = "nearestCancerGeneUpstream1";
my $svs_nearestCodingGeneDownstream1   = "nearestCodingGeneDownstream1";
my $svs_nearestCancerGeneDownstream1   = "nearestCancerGeneDownstream1";
my $svs_nearestCodingGeneUpstream2     = "nearestCodingGeneUpstream2";
my $svs_nearestCancerGeneUpstream2     = "nearestCancerGeneUpstream2";
my $svs_nearestCodingGeneDownstream2   = "nearestCodingGeneDownstream2";
my $svs_nearestCancerGeneDownstream2   = "nearestCancerGeneDownstream2";
my $svs_affectedCodingGenesTADestimate = "affectedCodingGenesTADestimate";
my $svs_svtype                         = "svtype";  
my $svs_eventInversion                 = "eventInversion";


foreach my $pid (@pids){

  $svs{$pid};
  my $sv_file = $mutation_files{$pid}{"sv_file"};
  my @svs_loci;
  my %header_index=();

  my $header_line=`head -1 $sv_file`;
  my @header = split ("\t", $header_line);
  my $header_elem_counter=0;
  foreach (@header){
    $header_index{$_}=$header_elem_counter;
    $header_elem_counter++;
  }

  @svs_loci = `grep -v "\#" $sv_file`;

  foreach my $locuss (@svs_loci){

    my @locus = split ("\t", $locuss);

    $sample_info{"SV - total"}{$pid} = $sample_info{"SV - total"}{$pid} + 1;
    my $sv_type= "$locus[$header_index{$svs_svtype}] - $locus[$header_index{$svs_eventInversion}]";
    $sample_info{"SV - $sv_type"}{$pid} = $sample_info{"SV - $sv_type"}{$pid} + 1;

    # direct 
    for my $direct_index($header_index{$svs_gene1},$header_index{$svs_gene2}){
      my @genes1_direct = split (",", $locus[$direct_index]);
      foreach my $gene1_direct (@genes1_direct){
        my @gene_split = split(/\|/,$gene1_direct); 
        $output_hash{$gene_split[0]}{"SV_direct"}{"."}{$pid} = $output_hash{$gene_split[0]}{"SV_direct"}{"."}{$pid} + 1 unless ($gene_split[0] eq ".");
        $output_hash{$gene_split[0]}{"SV_direct"}{$sv_type}{$pid} = $output_hash{$gene_split[0]}{"SV_direct"}{$sv_type}{$pid} + 1 unless ($gene_split[0] eq ".");
      }
    }

    # upstream/downstream 
    for my $close_index ($header_index{$svs_nearestCodingGeneUpstream1},$header_index{$svs_nearestCancerGeneUpstream1},$header_index{$svs_nearestCodingGeneDownstream1},$header_index{$svs_nearestCancerGeneDownstream1},$header_index{$svs_nearestCodingGeneUpstream2},$header_index{$svs_nearestCancerGeneUpstream2},$header_index{$svs_nearestCodingGeneDownstream2},$header_index{$svs_nearestCancerGeneDownstream2}){
      if (abs($locus[$close_index+1])<100000){
        $output_hash{$locus[$close_index]}{"SV_near"}{"."}{$pid} = $output_hash{$locus[$close_index]}{"SV_near"}{"."}{$pid} + 1 unless ($locus[$close_index] eq ".");
        $output_hash{$locus[$close_index]}{"SV_near"}{$sv_type}{$pid} = $output_hash{$locus[$close_index]}{"SV_near"}{$sv_type}{$pid} + 1 unless ($locus[$close_index] eq ".");
      }
    }

    # TAD 
    my @genes_tad = split (",", $locus[$header_index{$svs_affectedCodingGenesTADestimate}]);
    foreach my $tad_gene (@genes_tad){
      $output_hash{$tad_gene}{"SV_tad"}{"."}{$pid} = $output_hash{$tad_gene}{"SV_tad"}{"."}{$pid} + 1 unless ($tad_gene eq ".");
      $output_hash{$tad_gene}{"SV_tad"}{$sv_type}{$pid} = $output_hash{$tad_gene}{"SV_tad"}{$sv_type}{$pid} + 1 unless ($tad_gene eq ".");
    }
  }
}

warn "\# DONE!\n";
warn "\n\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n\n";

################
## PARSE SNVS ##
################

warn "\# PARSING: SNV files ...\n";

my $snvs_ANNOVAR_FUNCTION = "ANNOVAR_FUNCTION";
my $snvs_GENE = "GENE";
my $snvs_EXONIC_CLASSIFICATION = "EXONIC_CLASSIFICATION";

foreach my $pid (@pids){

  my $mutation_file = $mutation_files{$pid}{"mutation_file"};
  my $mutation_file_germline = $mutation_files{$pid}{"mutation_file_germline"};

  my @snvs_loci;
  my %snvs_header_index=();

  my $snvs_header_line=`head -1000 $mutation_file | grep "^\#CHROM"`;
  my @snvs_header = split ("\t", $snvs_header_line);
  my $snvs_header_elem_counter=0;
  foreach (@snvs_header){
    $snvs_header_index{$_}=$snvs_header_elem_counter;
    $snvs_header_elem_counter++;
  }

  # technically, the headers should be the same for somatic and germline files, but just to be sure:
  my @snvs_loci_germline;
  my %snvs_header_index_germline=();

  if (-e $mutation_file_germline){
    my $snvs_header_line_germline=`head -1000 $mutation_file_germline | grep "^\#CHROM"`;
    my @snvs_header_germline = split ("\t", $snvs_header_line_germline);
    my $snvs_header_elem_counter_germline=0;
    foreach (@snvs_header_germline){
      $snvs_header_index_germline{$_}=$snvs_header_elem_counter_germline;
      $snvs_header_elem_counter_germline++;
    }
  }

  @snvs_loci = `grep -v "^\#" $mutation_file | grep -vw intergenic`;
  
  foreach (@snvs_loci){
    my @snv_locus = split ("\t", $_);
    my $snv_string = $pid."\t".$snv_locus[$snvs_header_index{$snvs_ANNOVAR_FUNCTION}]."\t".$snv_locus[$snvs_header_index{$snvs_GENE}]."\t".$snv_locus[$snvs_header_index{$snvs_EXONIC_CLASSIFICATION}];

    chomp $snv_string;
    my @snv = split ("\t", $snv_string);
    my $pid = $snv[0];  
    # first check/split for ";"

    if (($snv[1] =~ tr/\;//) > 0){
      my @snv1 = split ("\;",$snv[1]);
      my @snv2 = split ("\;",$snv[2]);
      my @snv3 = split ("\;",$snv[3]);
      my $j = 0;
      for (my $i=0; $i <= length(@snv1); $i++){
        my @genes = split (",", $snv2[$i]);
        foreach my $gene (@genes){
          next if $gene =~ m/^ENST00/; # skip if reporting multiple transcripts per gene
          next if $snv1[$i] eq "intergenic";
          if (($gene =~ tr/\(//) > 0){
            $gene =~ m/^(.*?)\(.*$/;
            $gene = $1;
          } 

          # inititate_hash
          setup_gene_snv($gene, $pid) if (!(defined $output_hash{$gene}));

          if ($snv[1] eq "exonic"){
              $output_hash{$gene}{$snv[1]}{$snv[3]}{$pid} = $output_hash{$gene}{$snv[1]}{$snv[3]}{$pid} + 1;
          }
          else {
              $output_hash{$gene}{$snv[1]}{"."}{$pid} = $output_hash{$gene}{$snv[1]}{"."}{$pid} + 1;
          }
        }
        $j = $j + 1 if ($snv1[$i] eq "exonic" && length(@snv3)>1);
      }
    }
  
    # if no ";" then check for ",";
    else {
      next if $snv[1] eq "intergenic";
      my @genes = split (",", $snv[2]);
      foreach my $gene (@genes){
        next if $gene =~ m/^ENST00/; # skip if reporting multiple transcripts per gene
        if (($gene =~ tr/\(//) > 0){
          $gene =~ m/^(.*?)\(.*$/;
          $gene = $1;
        }
        # inititate_hash
        setup_gene_snv($gene, $pid) if (!(defined $output_hash{$gene}));

        $snv[3] = "." unless ( $snv[1] eq "exonic" );
        if ($snv[1] eq "exonic"){
            $output_hash{$gene}{$snv[1]}{$snv[3]}{$pid} = $output_hash{$gene}{$snv[1]}{$snv[3]}{$pid} + 1;
        }
        else {
            $output_hash{$gene}{$snv[1]}{"."}{$pid} = $output_hash{$gene}{$snv[1]}{"."}{$pid} + 1;
        }
      }
    }
  }

  if (-e $mutation_file_germline){
    @snvs_loci_germline = `grep -v "^\#" $mutation_file_germline | grep -vw intergenic`;
  }

  foreach (@snvs_loci_germline){
    my @snv_locus_germline = split ("\t", $_);
    my $snv_string_germline = $pid."\t".$snv_locus_germline[$snvs_header_index_germline{$snvs_ANNOVAR_FUNCTION}]."\t".$snv_locus_germline[$snvs_header_index_germline{$snvs_GENE}]."\t".$snv_locus_germline[$snvs_header_index_germline{$snvs_EXONIC_CLASSIFICATION}];

    chomp $snv_string_germline;
    my @snv_germline = split ("\t", $snv_string_germline);
    my $pid = $snv_germline[0];  
    # first check/split for ";"

    if (($snv_germline[1] =~ tr/\;//) > 0){
      my @snv1_germline = split ("\;",$snv_germline[1]);
      my @snv2_germline = split ("\;",$snv_germline[2]);
      my @snv3_germline = split ("\;",$snv_germline[3]);
      my $j = 0;
      for (my $i=0; $i <= length(@snv1_germline); $i++){
        my @genes_germline = split (",", $snv2_germline[$i]);
        foreach my $gene (@genes_germline){
          next if $gene =~ m/^ENST00/; # skip if reporting multiple transcripts per gene
          next if $snv1_germline[$i] eq "intergenic";
          if (($gene =~ tr/\(//) > 0){
            $gene =~ m/^(.*?)\(.*$/;
            $gene = $1;
          } 

          # inititate_hash
          setup_gene_snv_germline($gene, $pid) if (!(defined $output_hash_germline{$gene}));

          if ($snv_germline[1] eq "exonic"){
              $output_hash_germline{$gene}{$snv_germline[1]}{$snv_germline[3]}{$pid} = $output_hash_germline{$gene}{$snv_germline[1]}{$snv_germline[3]}{$pid} + 1;
          }
          else {
              $output_hash_germline{$gene}{$snv_germline[1]}{"."}{$pid} = $output_hash_germline{$gene}{$snv_germline[1]}{"."}{$pid} + 1;
          }
        }
        $j = $j + 1 if ($snv1_germline[$i] eq "exonic" && length(@snv3_germline)>1);
      }
    }
  
    # if no ";" then check for ",";
    else {
      next if $snv_germline[1] eq "intergenic";
      my @genes_germline = split (",", $snv_germline[2]);
      foreach my $gene (@genes_germline){
        next if $gene =~ m/^ENST00/; # skip if reporting multiple transcripts per gene
        if (($gene =~ tr/\(//) > 0){
          $gene =~ m/^(.*?)\(.*$/;
          $gene = $1;
        }
        # inititate_hash
        setup_gene_snv_germline($gene, $pid) if (!(defined $output_hash_germline{$gene}));

        $snv_germline[3] = "." unless ( $snv_germline[1] eq "exonic" );
        if ($snv_germline[1] eq "exonic"){
            $output_hash_germline{$gene}{$snv_germline[1]}{$snv_germline[3]}{$pid} = $output_hash_germline{$gene}{$snv_germline[1]}{$snv_germline[3]}{$pid} + 1;

        }
        else {
            $output_hash_germline{$gene}{$snv_germline[1]}{"."}{$pid} = $output_hash_germline{$gene}{$snv_germline[1]}{"."}{$pid} + 1;
        }
      }
    }
  }
}

warn "\# DONE!\n";
warn "\n\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n\n";

##################
## PARSE INDELs ##
##################

warn "\# PARSING: indel files ...\n";

my $indel_ANNOVAR_FUNCTION = "ANNOVAR_FUNCTION";
my $indel_GENE = "GENE";
my $indel_EXONIC_CLASSIFICATION = "EXONIC_CLASSIFICATION";

foreach my $pid (@pids){

  my $indel_file = $mutation_files{$pid}{"indel_file"};
  my $indel_file_germline = $mutation_files{$pid}{"indel_file_germline"};

  my @indels_loci;
  my %indels_header_index=();

  my $indels_header_line=`head -1000 $indel_file | grep "^\#CHROM"`;
  my @indels_header = split ("\t", $indels_header_line);
  my $indels_header_elem_counter=0;
  foreach (@indels_header){
    $indels_header_index{$_}=$indels_header_elem_counter;
    $indels_header_elem_counter++;
  }
  
  # technically, the headers should be the same for somatic and germline files, but just to be sure:
  my @indels_loci_germline;
  my %indels_header_index_germline=();

  if (-e $indel_file_germline){
    my $indels_header_line_germline=`head -1000 $indel_file_germline | grep "^\#CHROM"`;
    my @indels_header_germline = split ("\t", $indels_header_line_germline);
    my $indels_header_elem_counter_germline=0;
    foreach (@indels_header_germline){
      $indels_header_index_germline{$_}=$indels_header_elem_counter_germline;
      $indels_header_elem_counter_germline++;
    }
  }

  @indels_loci = `grep -v "^\#" $indel_file | grep -vw intergenic`;

  foreach (@indels_loci){
    my @indel_locus = split ("\t", $_);
    my $snv_string = $pid."\t".$indel_locus[$indels_header_index{$indel_ANNOVAR_FUNCTION}]."\t".$indel_locus[$indels_header_index{$indel_GENE}]."\t".$indel_locus[$indels_header_index{$indel_EXONIC_CLASSIFICATION}];

    chomp $snv_string;
    my @snv = split ("\t", $snv_string);
    my $pid = $snv[0];  

    # first check/split for ";"
    if (($snv[1] =~ tr/\;//) > 0){
      my @snv1 = split ("\;",$snv[1]);
      my @snv2 = split ("\;",$snv[2]);
      my @snv3 = split ("\;",$snv[3]);
      my $j = 0;
      for (my $i=0; $i <= length(@snv1); $i++){
        my @genes = split (",", $snv2[$i]);
        foreach my $gene (@genes){
          next if $gene =~ m/^ENST00/; # skip if reporting multiple transcripts per gene
          next if $snv1[$i] eq "intergenic";
          if (($gene =~ tr/\(//) > 0){
            $gene =~ m/^(.*?)\(.*$/;
            $gene = $1;
          } 

          # inititate_hash
          setup_gene_indel($gene, $pid) if (!(defined $output_hash{$gene}));
          if ($snv[1] eq "exonic"){
              $output_hash{$gene}{"INDEL_".$snv[1]}{$snv[3]}{$pid} = $output_hash{$gene}{"INDEL_".$snv[1]}{$snv[3]}{$pid} + 1;
          }
          else {
              $output_hash{$gene}{"INDEL_".$snv[1]}{"."}{$pid} = $output_hash{$gene}{"INDEL_".$snv[1]}{"."}{$pid} + 1;
          }
        }
        $j = $j + 1 if ($snv1[$i] eq "exonic" && length(@snv3)>1);
      }
    }
  
    # if no ";" then check for ",";
    else {
      next if $snv[1] eq "intergenic";
      my @genes = split (",", $snv[2]);
      foreach my $gene (@genes){
        next if $gene =~ m/^ENST00/; # skip if reporting multiple transcripts per gene
        if (($gene =~ tr/\(//) > 0){
          $gene =~ m/^(.*?)\(.*$/;
          $gene = $1;
        }
        # inititate_hash
        setup_gene_indel($gene, $pid) if (!(defined $output_hash{$gene}));
        if ($snv[1] eq "exonic"){
            $output_hash{$gene}{"INDEL_".$snv[1]}{$snv[3]}{$pid} = $output_hash{$gene}{"INDEL_".$snv[1]}{$snv[3]}{$pid} + 1;
        }
        else {
            $output_hash{$gene}{"INDEL_".$snv[1]}{"."}{$pid} = $output_hash{$gene}{"INDEL_".$snv[1]}{"."}{$pid} + 1;
        }
      }
    }
  }

  if (-e $indel_file_germline){
    @indels_loci_germline = `grep -v "^\#" $indel_file_germline | grep -vw intergenic`;
  }

  foreach (@indels_loci_germline){
    my @indel_locus_germline = split ("\t", $_);
    my $snv_string_germline = $pid."\t".$indel_locus_germline[$indels_header_index_germline{$indel_ANNOVAR_FUNCTION}]."\t".$indel_locus_germline[$indels_header_index_germline{$indel_GENE}]."\t".$indel_locus_germline[$indels_header_index_germline{$indel_EXONIC_CLASSIFICATION}];

    chomp $snv_string_germline;
    my @snv_germline = split ("\t", $snv_string_germline);
    my $pid = $snv_germline[0];  

    # first check/split for ";"
    if (($snv_germline[1] =~ tr/\;//) > 0){
      my @snv1_germline = split ("\;",$snv_germline[1]);
      my @snv2_germline = split ("\;",$snv_germline[2]);
      my @snv3_germline = split ("\;",$snv_germline[3]);
      my $j = 0;
      for (my $i=0; $i <= length(@snv1_germline); $i++){
        my @genes_germline = split (",", $snv2_germline[$i]);
        foreach my $gene (@genes_germline){
          next if $gene =~ m/^ENST00/; # skip if reporting multiple transcripts per gene
          next if $snv1_germline[$i] eq "intergenic";
          if (($gene =~ tr/\(//) > 0){
            $gene =~ m/^(.*?)\(.*$/;
            $gene = $1;
          } 

          # inititate_hash
          setup_gene_indel_germline($gene, $pid) if (!(defined $output_hash_germline{$gene}));
          if ($snv_germline[1] eq "exonic"){
              $output_hash_germline{$gene}{"INDEL_germline".$snv_germline[1]}{$snv_germline[3]}{$pid} = $output_hash_germline{$gene}{"INDEL_germline".$snv_germline[1]}{$snv_germline[3]}{$pid} + 1;
          }
          else {
              $output_hash_germline{$gene}{"INDEL_germline".$snv_germline[1]}{"."}{$pid} = $output_hash_germline{$gene}{"INDEL_germline".$snv_germline[1]}{"."}{$pid} + 1;
          }
        }
        $j = $j + 1 if ($snv1_germline[$i] eq "exonic" && length(@snv3_germline)>1);
      }
    }
  
    # if no ";" then check for ",";
    else {
      next if $snv_germline[1] eq "intergenic";
      my @genes_germline = split (",", $snv_germline[2]);
      foreach my $gene (@genes_germline){
        next if $gene =~ m/^ENST00/; # skip if reporting multiple transcripts per gene
        if (($gene =~ tr/\(//) > 0){
          $gene =~ m/^(.*?)\(.*$/;
          $gene = $1;
        }
        # inititate_hash
        setup_gene_indel_germline($gene, $pid) if (!(defined $output_hash_germline{$gene}));
        if ($snv_germline[1] eq "exonic"){
            $output_hash_germline{$gene}{"INDEL_germline_".$snv_germline[1]}{$snv_germline[3]}{$pid} = $output_hash_germline{$gene}{"INDEL_germline_".$snv_germline[1]}{$snv_germline[3]}{$pid} + 1;
        }
        else {
            $output_hash_germline{$gene}{"INDEL_germline_".$snv_germline[1]}{"."}{$pid} = $output_hash_germline{$gene}{"INDEL_germline_".$snv_germline[1]}{"."}{$pid} + 1;
        }
      }
    }
  }
}

warn "\# DONE!\n";
warn "\n\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n\n";

####################
## Fusion genes
####################
# Fusion genes for available RNAseq samples
## 

warn "\# PARSING FUSION GENES.... \n";
for my $pid (@pids) {

  if($mutation_files{$pid}{'rnaseq_fusion_file'} ne '') {
    open(IN, $mutation_files{$pid}{'rnaseq_fusion_file'}) || die "cant open the file fusion file\n";

    while(<IN>) {
      chomp;
      if($_!~/^#gene1/ && $_=~/splice-site/ && $_!~/read-through/) {
        my @ss = split(/\t/, $_);
        if($ss[14] !~ /low/) {
          my @split_gene_1 = split(/,/, $ss[0]);
          my @split_gene_2 = split(/,/, $ss[1]);
          push(@split_gene_1, @split_gene_2);

          my @split_gene_clean = map{my ($foo)=$_=~s/\(\d+\)//; $foo} @split_gene_1;

          for my $gene (@split_gene_1) {
            if(defined $output_hash{$gene}{"rna_fusion"}{"."}{$pid}) {
              $output_hash{$gene}{"rna_fusion"}{"."}{$pid} += 1;
            }
            else {
              $output_hash{$gene}{"rna_fusion"}{"."}{$pid} = 1;
            }
          }
        }
      }
    }
  }
}

warn "\# DONE!\n";
warn "\n\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n\n";

####################
## OPEN OUT FILES ##
####################

open (my $outfile_fh, ">$outfile")                         or die "\# ERROR: cannot open file $outfile\n\n";
open (my $outfile_dupl_fh, ">$outfile_dupl")               or die "\# ERROR: cannot open file $outfile_dupl\n\n";
open (my $outfile_sample_info_fh, ">$outfile_sample_info") or die "\# ERROR: cannot open file $outfile_sample_info\n\n";
open (my $outfile_kataegis_fh, ">$outfile_kataegis")       or die "\# ERROR: cannot open file $outfile_kataegis\n\n";

#######################
## PRINT SAMPLE_INFO ##
#######################

warn "\# PREPARING SAMPLE INFO TABLE ...\n";

foreach my $pid (@pids){
  chomp $pid;
  print $outfile_sample_info_fh "\t$pid";
}
print $outfile_sample_info_fh "\n";

foreach my $feature (sort keys %sample_info) {
  print $outfile_sample_info_fh  "$feature";
  foreach my $pid (@pids){
    print $outfile_sample_info_fh "\t$sample_info{$feature}{$pid}";
  }
  print $outfile_sample_info_fh "\n";
}

warn "\# DONE!\n";
warn "\n\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n\n";

###############
## PRINT OUT ##
###############

# put genes with duplicated gene names in separate list:
my %duplication_list=();

warn "\# PREPARING OUTPUT TABLE ...\n";

foreach my $pid (@pids){
  chomp $pid;
  print $outfile_fh "\t$pid";
}
print $outfile_fh "\n";

foreach my $pid (@pids){
  chomp $pid;
  print $outfile_dupl_fh "\t$pid";
}
print $outfile_dupl_fh "\n";

foreach my $gene (sort keys %output_hash) {
  if ($genes{$gene}{"count"} > 1){
    $duplication_list{$gene}=1;
    warn "\nWARNING: \"$gene\" will be in separate table, as the gene name is defined multiple times in the gene models\n";
    #next;
  }

  if (defined $exclusion_list{$gene}){
    warn "\nWARNING: skipping \"$gene\" as it is defined in the exclude list\n";
    next;
  }
#  next if (defined $exclusion_list{$gene});
    
  my $outline = "";
  my $count = 0; 
 
  $outline = $outline."$gene\t";

  # precount

  my %pid_precount;
  my %pid_precount_germline;
  my %pid_precount_svs;
  my %pid_precount_cnvs;
  my %pid_precount_fusion;

  foreach my $pid (@pids){ 
    $pid_precount{$pid} = 1 if ($output_hash{$gene}{"exonic"}{"nonsynonymous SNV"}{$pid}) ge 1;
    $pid_precount{$pid} = 1 if ($output_hash{$gene}{"exonic"}{"stoploss"}{$pid}) ge 1;
    $pid_precount{$pid} = 1 if ($output_hash{$gene}{"exonic"}{"stopgain"}{$pid}) ge 1;
    $pid_precount{$pid} = 1 if ($output_hash{$gene}{"splicing"}{"."}{$pid}) ge 1;
    $pid_precount{$pid} = 1 if ($output_hash{$gene}{"INDEL_exonic"}{"frameshift deletion"}{$pid}) ge 1;
    $pid_precount{$pid} = 1 if ($output_hash{$gene}{"INDEL_exonic"}{"frameshift insertion"}{$pid}) ge 1;
    $pid_precount{$pid} = 1 if ($output_hash{$gene}{"INDEL_exonic"}{"nonframeshift deletion"}{$pid}) ge 1;
    $pid_precount{$pid} = 1 if ($output_hash{$gene}{"INDEL_exonic"}{"nonframeshift insertion"}{$pid}) ge 1;
    $pid_precount{$pid} = 1 if ($output_hash{$gene}{"INDEL_exonic"}{"stopgain"}{$pid}) ge 1;
    $pid_precount{$pid} = 1 if ($output_hash{$gene}{"INDEL_exonic"}{"stoploss"}{$pid}) ge 1;
    $pid_precount{$pid} = 1 if ($output_hash{$gene}{"INDEL_splicing"}{"."}{$pid}) ge 1;
    $pid_precount{$pid} = 1 if ($output_hash{$gene}{"INDEL_ncRNA_exonic"}{"."}{$pid}) ge 1;
    $pid_precount{$pid} = 1 if ($output_hash{$gene}{"ncRNA_exonic"}{"."}{$pid}) ge 1;
    $pid_precount{$pid} = 1 if ($output_hash{$gene}{"HOMODEL"}{"."}{$pid}) ge 1;
    $pid_precount{$pid} = 1 if ($output_hash{$gene}{"GAIN"}{"."}{$pid}) ge 1;
    $pid_precount{$pid} = 1 if ($output_hash{$gene}{"HIGHGAIN"}{"."}{$pid}) ge 1;
    $pid_precount{$pid} = 1 if ($output_hash{$gene}{"LOSS"}{"."}{$pid}) ge 1;
    $pid_precount{$pid} = 1 if ($output_hash{$gene}{"LOH"}{"."}{$pid}) ge 1;
    $pid_precount{$pid}  = 1 if ($output_hash{$gene}{"SV_direct"}{"."}{$pid}) ge 1;
    $pid_precount_svs{$pid}  = 1 if ($pid_precount{$pid} ge 1); 
    $pid_precount_svs{$pid}  = 1 if ($output_hash{$gene}{"SV_direct"}{"."}{$pid}) ge 1;
    $pid_precount_cnvs{$pid} = 1 if ($output_hash{$gene}{"HOMODEL"}{"."}{$pid}) ge 1;
    $pid_precount_cnvs{$pid} = 1 if ($output_hash{$gene}{"GAIN"}{"."}{$pid}) ge 1;
    $pid_precount_cnvs{$pid} = 1 if ($output_hash{$gene}{"HIGHGAIN"}{"."}{$pid}) ge 1;
    $pid_precount_cnvs{$pid} = 1 if ($output_hash{$gene}{"LOSS"}{"."}{$pid}) ge 1;
    $pid_precount_cnvs{$pid} = 1 if ($output_hash{$gene}{"LOH"}{"."}{$pid}) ge 1;
    $pid_precount_fusion{$pid} = 1 if ($output_hash{$gene}{"rna_fusion"}{"."}{$pid}) ge 1;
  }

  my $pids_hit = scalar(keys %pid_precount);
  my $pids_hit_svs = scalar(keys %pid_precount_svs);
  my $pids_hit_cnv = scalar(keys %pid_precount_cnvs);

  foreach my $pid (@pids){

    my $pidline = "";
    my $mut_SYN= ($output_hash{$gene}{"exonic"}{"synonymous SNV"}{$pid}) ge 1;
    my $mut_SYN_germline= ($output_hash_germline{$gene}{"exonic"}{"synonymous SNV"}{$pid}) ge 1;
    my $mut_NONSYN= ($output_hash{$gene}{"exonic"}{"nonsynonymous SNV"}{$pid}) ge 1;
    my $mut_NONSYN_germline= ($output_hash_germline{$gene}{"exonic"}{"nonsynonymous SNV"}{$pid}) ge 1;
    my $mut_STOPL = ($output_hash{$gene}{"exonic"}{"stoploss"}{$pid}) ge 1;
    my $mut_STOPL_germline = ($output_hash_germline{$gene}{"exonic"}{"stoploss"}{$pid}) ge 1;
    my $mut_STOPG = ($output_hash{$gene}{"exonic"}{"stopgain"}{$pid}) ge 1;
    my $mut_STOPG_germline = ($output_hash_germline{$gene}{"exonic"}{"stopgain"}{$pid}) ge 1;
    my $mut_SPLICE= ($output_hash{$gene}{"splicing"}{"."}{$pid}) ge 1;
    my $mut_SPLICE_germline= ($output_hash_germline{$gene}{"splicing"}{"."}{$pid}) ge 1;
    my $mut_5UTR  = ($output_hash{$gene}{"UTR5"}{"."}{$pid}) ge 1;
    my $mut_5UTR_germline  = ($output_hash_germline{$gene}{"UTR5"}{"."}{$pid}) ge 1;
    my $mut_UP    = ($output_hash{$gene}{"upstream"}{"."}{$pid}) ge 1;
    my $mut_UP_germline    = ($output_hash_germline{$gene}{"upstream"}{"."}{$pid}) ge 1;
    my $mut_3UTR  = ($output_hash{$gene}{"UTR3"}{"."}{$pid}) ge 1;
    my $mut_3UTR_germline  = ($output_hash_germline{$gene}{"UTR3"}{"."}{$pid}) ge 1;
    my $mut_DOWN  = ($output_hash{$gene}{"downstream"}{"."}{$pid}) ge 1;
    my $mut_DOWN_germline  = ($output_hash_germline{$gene}{"downstream"}{"."}{$pid}) ge 1;
    my $mut_NCRNA = ($output_hash{$gene}{"ncRNA_exonic"}{"."}{$pid}) ge 1;
    my $mut_NCRNA_germline = ($output_hash_germline{$gene}{"ncRNA_exonic"}{"."}{$pid}) ge 1;
    my $mut_INTRON= ($output_hash{$gene}{"intronic"}{"."}{$pid}) ge 1;
    my $mut_INTRON_germline= ($output_hash_germline{$gene}{"intronic"}{"."}{$pid}) ge 1;

    my $mut_SHM   = ($output_hash{$gene}{"SHM"}{"."}{$pid}) ge 1;

    my $sv_DIRECT = ($output_hash{$gene}{"SV_direct"}{"."}{$pid}) ge 1;
    my $sv_NEAR   = ($output_hash{$gene}{"SV_near"}{"."}{$pid}) ge 1;
    my $sv_TAD    = ($output_hash{$gene}{"SV_tad"}{"."}{$pid}) ge 1;

    my $cnv_HIGHGAIN=($output_hash{$gene}{"HIGHGAIN"}{"."}{$pid}) ge 1;
    my $cnv_GAIN   = ($output_hash{$gene}{"GAIN"}{"."}{$pid}) ge 1;
    my $cnv_LOSS   = ($output_hash{$gene}{"LOSS"}{"."}{$pid}) ge 1;
    my $cnv_HOMODEL= ($output_hash{$gene}{"HOMODEL"}{"."}{$pid}) ge 1;
    my $cnv_LOH    = ($output_hash{$gene}{"LOH"}{"."}{$pid})  ge 1;

    my $cnv_chrGAIN   = ($output_hash{$gene}{"chrGAIN"}{"."}{$pid}) ge 1;
    my $cnv_chrLOSS   = ($output_hash{$gene}{"chrLOSS"}{"."}{$pid}) ge 1;
    my $cnv_chrHOMODEL= ($output_hash{$gene}{"chrHOMODEL"}{"."}{$pid}) ge 1;
    my $cnv_chrLOH    = ($output_hash{$gene}{"chrLOH"}{"."}{$pid})  ge 1;

    my $indel_FSD   = ($output_hash{$gene}{"INDEL_exonic"}{"frameshift deletion"}{$pid}) ge 1;
    my $indel_FSD_germline   = ($output_hash_germline{$gene}{"INDEL_exonic"}{"frameshift deletion"}{$pid}) ge 1;
    my $indel_FSI   = ($output_hash{$gene}{"INDEL_exonic"}{"frameshift insertion"}{$pid}) ge 1;
    my $indel_FSI_germline   = ($output_hash_germline{$gene}{"INDEL_exonic"}{"frameshift insertion"}{$pid}) ge 1;
    my $indel_NFSD  = ($output_hash{$gene}{"INDEL_exonic"}{"nonframeshift deletion"}{$pid}) ge 1;
    my $indel_NFSD_germline  = ($output_hash_germline{$gene}{"INDEL_exonic"}{"nonframeshift deletion"}{$pid}) ge 1;
    my $indel_NFSI  = ($output_hash{$gene}{"INDEL_exonic"}{"nonframeshift insertion"}{$pid}) ge 1;
    my $indel_NFSI_germline  = ($output_hash_germline{$gene}{"INDEL_exonic"}{"nonframeshift insertion"}{$pid}) ge 1;
    my $indel_STOPG = ($output_hash{$gene}{"INDEL_exonic"}{"stopgain"}{$pid}) ge 1;
    my $indel_STOPG_germline = ($output_hash_germline{$gene}{"INDEL_exonic"}{"stopgain"}{$pid}) ge 1;
    my $indel_STOPL = ($output_hash{$gene}{"INDEL_exonic"}{"stoploss"}{$pid}) ge 1;
    my $indel_STOPL_germline = ($output_hash_germline{$gene}{"INDEL_exonic"}{"stoploss"}{$pid}) ge 1;
    my $indel_SPLICE= ($output_hash{$gene}{"INDEL_splicing"}{"."}{$pid}) ge 1;
    my $indel_SPLICE_germline= ($output_hash_germline{$gene}{"INDEL_splicing"}{"."}{$pid}) ge 1;
    my $indel_NCRNA = ($output_hash{$gene}{"INDEL_ncRNA_exonic"}{"."}{$pid}) ge 1;
    my $indel_NCRNA_germline = ($output_hash_germline{$gene}{"INDEL_ncRNA_exonic"}{"."}{$pid}) ge 1;
    my $indel_UP    = ($output_hash{$gene}{"INDEL_upstream"}{"."}{$pid}) ge 1;
    my $indel_UP_germline    = ($output_hash_germline{$gene}{"INDEL_upstream"}{"."}{$pid}) ge 1;
    my $indel_5UTR  = ($output_hash{$gene}{"INDEL_UTR5"}{"."}{$pid}) ge 1;
    my $indel_5UTR_germline  = ($output_hash_germline{$gene}{"INDEL_UTR5"}{"."}{$pid}) ge 1;
    my $indel_DOWN  = ($output_hash{$gene}{"INDEL_downstream"}{"."}{$pid}) ge 1;
    my $indel_DOWN_germline  = ($output_hash_germline{$gene}{"INDEL_downstream"}{"."}{$pid}) ge 1;
    my $indel_3UTR  = ($output_hash{$gene}{"INDEL_UTR3"}{"."}{$pid}) ge 1;
    my $indel_3UTR_germline  = ($output_hash_germline{$gene}{"INDEL_UTR3"}{"."}{$pid}) ge 1;
    my $indel_INTRON= ($output_hash{$gene}{"INDEL_intronic"}{"."}{$pid}) ge 1;
    my $indel_INTRON_germline= ($output_hash_germline{$gene}{"INDEL_intronic"}{"."}{$pid}) ge 1;

    my $rna_fusion = ($output_hash{$gene}{"rna_fusion"}{"."}{$pid}) ge 1; 

      $pidline = $pidline. "synonymous_SNV:".$output_hash{$gene}{"exonic"}{"synonymous SNV"}{$pid}."\;" if $mut_SYN;
      $pidline = $pidline. "synonymous_SNV_germline:".$output_hash_germline{$gene}{"exonic"}{"synonymous SNV"}{$pid}."\;" if $mut_SYN_germline;
      $pidline = $pidline. "nonsynonymous_SNV:".$output_hash{$gene}{"exonic"}{"nonsynonymous SNV"}{$pid}."\;" if $mut_NONSYN;
      $pidline = $pidline. "nonsynonymous_SNV_germline:".$output_hash_germline{$gene}{"exonic"}{"nonsynonymous SNV"}{$pid}."\;" if $mut_NONSYN_germline;
      $pidline = $pidline. "stoploss_snv:".$output_hash{$gene}{"exonic"}{"stoploss"}{$pid}."\;" if ($mut_STOPL);
      $pidline = $pidline. "stoploss_snv_germline:".$output_hash_germline{$gene}{"exonic"}{"stoploss"}{$pid}."\;" if ($mut_STOPL_germline);
      $pidline = $pidline. "stopgain_snv:".$output_hash{$gene}{"exonic"}{"stopgain"}{$pid}."\;" if ($mut_STOPG);
      $pidline = $pidline. "stopgain_snv_germline:".$output_hash_germline{$gene}{"exonic"}{"stopgain"}{$pid}."\;" if ($mut_STOPG_germline);
      $pidline = $pidline. "splicing_snv:".$output_hash{$gene}{"splicing"}{"."}{$pid}."\;" if ($mut_SPLICE);
      $pidline = $pidline. "splicing_snv_germline:".$output_hash_germline{$gene}{"splicing"}{"."}{$pid}."\;" if ($mut_SPLICE_germline);

      $pidline = $pidline. "frameshift_deletion:".$output_hash{$gene}{"INDEL_exonic"}{"frameshift deletion"}{$pid}."\;"     if ($indel_FSD);
      $pidline = $pidline. "frameshift_deletion_germline:".$output_hash_germline{$gene}{"INDEL_exonic"}{"frameshift deletion"}{$pid}."\;"    if ($indel_FSI_germline);
      $pidline = $pidline. "frameshift_insertion:".$output_hash{$gene}{"INDEL_exonic"}{"frameshift insertion"}{$pid}."\;"    if ($indel_FSI);
      $pidline = $pidline. "frameshift_insertion_germline:".$output_hash_germline{$gene}{"INDEL_exonic"}{"frameshift insertion"}{$pid}."\;"    if ($indel_FSI_germline);
      $pidline = $pidline. "nonframeshift_deletion:".$output_hash{$gene}{"INDEL_exonic"}{"nonframeshift deletion"}{$pid}."\;"  if ($indel_NFSD);
      $pidline = $pidline. "nonframeshift_deletion_germline:".$output_hash_germline{$gene}{"INDEL_exonic"}{"nonframeshift deletion"}{$pid}."\;"  if ($indel_NFSD_germline);
      $pidline = $pidline. "nonframeshift_insertion:".$output_hash{$gene}{"INDEL_exonic"}{"nonframeshift insertion"}{$pid}."\;" if ($indel_NFSI);
      $pidline = $pidline. "nonframeshift_insertion_germline:".$output_hash_germline{$gene}{"INDEL_exonic"}{"nonframeshift insertion"}{$pid}."\;" if ($indel_NFSI_germline);

      $pidline = $pidline. "stoploss_indel:".$output_hash{$gene}{"INDEL_exonic"}{"stoploss"}{$pid}."\;" if ($indel_STOPL);
      $pidline = $pidline. "stoploss_indel_germline:".$output_hash_germline{$gene}{"INDEL_exonic"}{"stoploss"}{$pid}."\;" if ($indel_STOPL_germline);
      $pidline = $pidline. "stopgain_indel:".$output_hash{$gene}{"INDEL_exonic"}{"stopgain"}{$pid}."\;" if ($indel_STOPG);
      $pidline = $pidline. "stopgain_indel_germline:".$output_hash_germline{$gene}{"INDEL_exonic"}{"stopgain"}{$pid}."\;" if ($indel_STOPG_germline);
      $pidline = $pidline. "splicing_indel:".$output_hash{$gene}{"INDEL_exonic"}{"splicing"}{$pid}."\;" if ($indel_SPLICE);
      $pidline = $pidline. "splicing_indel_germline:".$output_hash_germline{$gene}{"INDEL_exonic"}{"splicing"}{$pid}."\;" if ($indel_SPLICE_germline);

      # always report direct non arm level CNVs
      $pidline = $pidline. "high_amplification:".$output_hash{$gene}{"HIGHGAIN"}{"."}{$pid}."\;" if ($cnv_HIGHGAIN && !$cnv_LOH); # changed from $pidline = $pidline. "high_amplification:".$output_hash{$gene}{"HIGHGAIN"}{"."}{$pid}."\;" if ($cnv_HIGHGAIN); 
      $pidline = $pidline. "highAmp_LOH:".$output_hash{$gene}{"LOH"}{"."}{$pid}."\;" if ($cnv_HIGHGAIN && $cnv_LOH); # added (see line above)
      $pidline = $pidline. "amplification:".$output_hash{$gene}{"GAIN"}{"."}{$pid}."\;" if ($cnv_GAIN && !$cnv_LOH); # changed from $pidline = $pidline. "amplification:".$output_hash{$gene}{"GAIN"}{"."}{$pid}."\;" if ($cnv_GAIN); 
      $pidline = $pidline. "amp_LOH:".$output_hash{$gene}{"LOH"}{"."}{$pid}."\;" if ($cnv_GAIN && $cnv_LOH); # added (see line above)
      $pidline = $pidline. "deletion:".$output_hash{$gene}{"LOSS"}{"."}{$pid}."\;" if ($cnv_LOSS && !$cnv_LOH); # changed from $pidline = $pidline. "deletion:".     $output_hash{$gene}{"LOSS"}{"."}{$pid}."\;" if ($cnv_LOSS); 
      $pidline = $pidline. "del_LOH:".$output_hash{$gene}{"LOH"}{"."}{$pid}."\;" if ($cnv_LOSS && $cnv_LOH); # added (see line above)
      $pidline = $pidline. "homo_del:".$output_hash{$gene}{"HOMODEL"}{"."}{$pid}."\;" if ($cnv_HOMODEL);
      $pidline = $pidline. "cn_LOH:".$output_hash{$gene}{"LOH"}{"."}{$pid}. "\;" if ($cnv_LOH && !$cnv_HIGHGAIN && !$cnv_GAIN && !$cnv_LOSS); # changed from $pidline = $pidline. "LOH:".        $output_hash{$gene}{"LOH"}{"."}{$pid}. "\;" if ($cnv_LOH);

      # report chr level event if also mutated??
      #$pidline = $pidline. "chrAmplification:".$output_hash{$gene}{"chrGAIN"}{"."}{$pid}."\;" if ($cnv_chrGAIN && !$cnv_chrLOH && ($indel_SPLICE||$indel_STOPG||$indel_STOPL||$indel_NFSI||$indel_NFSD||$indel_FSI||$indel_FSD||$mut_SPLICE||$mut_STOPG||$mut_STOPL||$mut_NONSYN||$sv_DIRECT));
      #$pidline = $pidline. "chrAmpLOH:".     $output_hash{$gene}{"chrLOH"}{"."}{$pid}."\;" if ($cnv_chrGAIN && $cnv_chrLOH && ($indel_SPLICE||$indel_STOPG||$indel_STOPL||$indel_NFSI||$indel_NFSD||$indel_FSI||$indel_FSD||$mut_SPLICE||$mut_STOPG||$mut_STOPL||$mut_NONSYN||$sv_DIRECT));
      #$pidline = $pidline. "chrDeletion:".     $output_hash{$gene}{"chrLOSS"}{"."}{$pid}."\;" if ($cnv_chrLOSS && !$cnv_chrLOH && ($indel_SPLICE||$indel_STOPG||$indel_STOPL||$indel_NFSI||$indel_NFSD||$indel_FSI||$indel_FSD||$mut_SPLICE||$mut_STOPG||$mut_STOPL||$mut_NONSYN||$sv_DIRECT));
      #$pidline = $pidline. "chrDeletionLOH:".     $output_hash{$gene}{"chrLOH"}{"."}{$pid}."\;" if ($cnv_chrLOSS && $cnv_chrLOH && ($indel_SPLICE||$indel_STOPG||$indel_STOPL||$indel_NFSI||$indel_NFSD||$indel_FSI||$indel_FSD||$mut_SPLICE||$mut_STOPG||$mut_STOPL||$mut_NONSYN||$sv_DIRECT));
      #$pidline = $pidline. "chrHomoDel:".      $output_hash{$gene}{"chrHOMODEL"}{"."}{$pid}."\;" if ($cnv_chrHOMODEL && ($indel_SPLICE||$indel_STOPG||$indel_STOPL||$indel_NFSI||$indel_NFSD||$indel_FSI||$indel_FSD||$mut_SPLICE||$mut_STOPG||$mut_STOPL||$mut_NONSYN||$sv_DIRECT));
      #$pidline = $pidline. "chrCnLOH:".        $output_hash{$gene}{"chrLOH"}{"."}{$pid}. "\;" if ($cnv_chrLOH && !$cnv_chrGAIN && !$cnv_chrLOSS && ($indel_SPLICE||$indel_STOPG||$indel_STOPL||$indel_NFSI||$indel_NFSD||$indel_FSI||$indel_FSD||$mut_SPLICE||$mut_STOPG||$mut_STOPL||$mut_NONSYN||$sv_DIRECT));

      $pidline = $pidline. "SV_direct:".$output_hash{$gene}{"SV_direct"}{"."}{$pid}."\;" if ($sv_DIRECT);
      $pidline = $pidline. "SV_near:"  .$output_hash{$gene}{"SV_near"}{"."}{$pid}  ."\;" if ($sv_NEAR && !($sv_DIRECT) &&                ($pids_hit >=1               || $pids_hit_svs >= $min_recurrence));
      $pidline = $pidline. "SV_TAD:"   .$output_hash{$gene}{"SV_tad"}{"."}{$pid}   ."\;" if ($sv_TAD  && !($sv_DIRECT) && !($sv_NEAR) && ($pids_hit >=$min_recurrence || $pids_hit_svs >= $min_recurrence));

      $pidline = $pidline. "ncRNA_exonic_snv:".$output_hash{$gene}{"ncRNA_exonic"}{"."}{$pid}."\;"         if ($mut_NCRNA);
      $pidline = $pidline. "ncRNA_exonic_snv_germline:".$output_hash_germline{$gene}{"ncRNA_exonic"}{"."}{$pid}."\;"         if ($mut_NCRNA_germline);
      $pidline = $pidline. "ncRNA_exonic_indel:".$output_hash{$gene}{"INDEL_ncRNA_exonic"}{"."}{$pid}."\;" if ($indel_NCRNA);
      $pidline = $pidline. "ncRNA_exonic_indel_germline:".$output_hash_germline{$gene}{"INDEL_ncRNA_exonic"}{"."}{$pid}."\;" if ($indel_NCRNA_germline);

      $pidline = $pidline. "UTR_5_snv:".  $output_hash{$gene}{"UTR5"}{"."}{$pid}."\;"       if $mut_5UTR;
      $pidline = $pidline. "UTR_5_snv_germline:".  $output_hash_germline{$gene}{"UTR5"}{"."}{$pid}."\;"       if $mut_5UTR_germline;
      $pidline = $pidline. "UTR_3_snv:".  $output_hash{$gene}{"UTR3"}{"."}{$pid}."\;"       if $mut_3UTR;
      $pidline = $pidline. "UTR_3_snv_germline:".  $output_hash_germline{$gene}{"UTR3"}{"."}{$pid}."\;"       if $mut_3UTR_germline;

      $pidline = $pidline. "UTR_5_indel:".$output_hash{$gene}{"INDEL_UTR5"}{"."}{$pid}."\;" if $indel_5UTR;
      $pidline = $pidline. "UTR_5_indel_germline:".$output_hash_germline{$gene}{"INDEL_UTR5"}{"."}{$pid}."\;" if $indel_5UTR_germline;
      $pidline = $pidline. "UTR_3_indel:".$output_hash{$gene}{"INDEL_UTR3"}{"."}{$pid}."\;" if $indel_3UTR;
      $pidline = $pidline. "UTR_3_indel_germline:".$output_hash_germline{$gene}{"INDEL_UTR3"}{"."}{$pid}."\;" if $indel_3UTR_germline;

      $pidline = $pidline. "intronic_snv:"  .$output_hash{$gene}{"intronic"}{"."}{$pid}."\;"       if $mut_INTRON;
      $pidline = $pidline. "intronic_snv_germline:"  .$output_hash_germline{$gene}{"intronic"}{"."}{$pid}."\;"       if $mut_INTRON_germline;
      $pidline = $pidline. "intronic_indel:".$output_hash{$gene}{"INDEL_intronic"}{"."}{$pid}."\;" if $indel_INTRON;
      $pidline = $pidline. "intronic_indel_germline:".$output_hash_germline{$gene}{"INDEL_intronic"}{"."}{$pid}."\;" if $indel_INTRON_germline;

      $pidline = $pidline. "upstream_snv:".$output_hash{$gene}{"upstream"}{"."}{$pid}."\;" if $mut_UP;
      $pidline = $pidline. "upstream_snv_germline:".$output_hash_germline{$gene}{"upstream"}{"."}{$pid}."\;" if $mut_UP_germline;
      $pidline = $pidline. "upstream_indel:".$output_hash{$gene}{"INDEL_upstream"}{"."}{$pid}."\;" if $indel_UP;
      $pidline = $pidline. "upstream_indel_germline:".$output_hash_germline{$gene}{"INDEL_upstream"}{"."}{$pid}."\;" if $indel_UP_germline;

      $pidline = $pidline. "rna_fusion:1;" if $rna_fusion;
      $pidline = $pidline. "kataegis:1;" if $mut_SHM;

    $outline = $outline.$pidline;
    $outline = $outline."\t" unless ($pids[(scalar @pids)-1] eq $pid);
    $count = $count + 1 if (length $pidline > 2);
  }
  if (defined $duplication_list{$gene}){
    print $outfile_dupl_fh "$outline\n" if ($count >= $min_recurrence); 
  }
  else {
    print $outfile_fh "$outline\n" if ($count >= $min_recurrence); 
  }
}

warn "\# DONE!\n";
warn "\n\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n\n";

close ($outfile_fh);
close ($outfile_dupl_fh);
close ($outfile_sample_info_fh);
close ($outfile_kataegis_fh);

exit;

#### SUBROUTINES #####

sub setup_gene_snv {
    my $gene = shift;
    my $pid = shift;
    $output_hash{$gene}{"exonic"}{"synonymous SNV"}{$pid} = 0;
    $output_hash{$gene}{"exonic"}{"nonsynonymous SNV"}{$pid} = 0;
    $output_hash{$gene}{"exonic"}{"stoploss"}{$pid} = 0;
    $output_hash{$gene}{"exonic"}{"stopgain"}{$pid} = 0;
    $output_hash{$gene}{"exonic"}{"unknown"}{$pid} = 0;
    $output_hash{$gene}{"intronic"}{"."}{$pid} = 0;
    $output_hash{$gene}{"ncRNA_exonic"}{"."}{$pid} = 0;
    $output_hash{$gene}{"ncRNA_intronic"}{"."}{$pid} = 0;
    $output_hash{$gene}{"ncRNA_splicing"}{"."}{$pid} = 0;
    $output_hash{$gene}{"splicing"}{"."}{$pid} = 0;
    $output_hash{$gene}{"upstream"}{"."}{$pid} = 0;
    $output_hash{$gene}{"downstream"}{"."}{$pid} = 0;
    $output_hash{$gene}{"UTR3"}{"."}{$pid} = 0;
    $output_hash{$gene}{"UTR5"}{"."}{$pid} = 0;
    $output_hash{$gene}{"SHM"}{"."}{$pid} = 0;
    $output_hash{$gene}{"GAIN"}{"."}{$pid} = 0;
    $output_hash{$gene}{"HIGHGAIN"}{"."}{$pid} = 0;
    $output_hash{$gene}{"LOSS"}{"."}{$pid} = 0;
    $output_hash{$gene}{"HOMODEL"}{"."}{$pid} = 0;
    $output_hash{$gene}{"LOH"}{"."}{$pid} = 0;
    $output_hash{$gene}{"chrGAIN"}{"."}{$pid} = 0;
    $output_hash{$gene}{"chrLOSS"}{"."}{$pid} = 0;
    $output_hash{$gene}{"chrHOMODEL"}{"."}{$pid} = 0;
    $output_hash{$gene}{"chrLOH"}{"."}{$pid} = 0;
}

sub setup_gene_snv_germline {
    my $gene = shift;
    my $pid = shift;
    $output_hash_germline{$gene}{"exonic"}{"synonymous SNV"}{$pid} = 0;
    $output_hash_germline{$gene}{"exonic"}{"nonsynonymous SNV"}{$pid} = 0;
    $output_hash_germline{$gene}{"exonic"}{"stoploss"}{$pid} = 0;
    $output_hash_germline{$gene}{"exonic"}{"stopgain"}{$pid} = 0;
    $output_hash_germline{$gene}{"exonic"}{"unknown"}{$pid} = 0;
    $output_hash_germline{$gene}{"intronic"}{"."}{$pid} = 0;
    $output_hash_germline{$gene}{"ncRNA_exonic"}{"."}{$pid} = 0;
    $output_hash_germline{$gene}{"ncRNA_intronic"}{"."}{$pid} = 0;
    $output_hash_germline{$gene}{"ncRNA_splicing"}{"."}{$pid} = 0;
    $output_hash_germline{$gene}{"splicing"}{"."}{$pid} = 0;
    $output_hash_germline{$gene}{"upstream"}{"."}{$pid} = 0;
    $output_hash_germline{$gene}{"downstream"}{"."}{$pid} = 0;
    $output_hash_germline{$gene}{"UTR3"}{"."}{$pid} = 0;
    $output_hash_germline{$gene}{"UTR5"}{"."}{$pid} = 0;
    $output_hash_germline{$gene}{"SHM"}{"."}{$pid} = 0;
}

sub setup_gene_indel {
    my $gene = shift;
    my $pid = shift;
    $output_hash{$gene}{"INDEL_exonic"}{"frameshift deletion"}{$pid} = 0;
    $output_hash{$gene}{"INDEL_exonic"}{"frameshift insertion"}{$pid} = 0;
    $output_hash{$gene}{"INDEL_exonic"}{"nonframeshift deletion"}{$pid} = 0;
    $output_hash{$gene}{"INDEL_exonic"}{"nonframeshift insertion"}{$pid} = 0;
    $output_hash{$gene}{"INDEL_exonic"}{"stopgain"}{$pid} = 0;
    $output_hash{$gene}{"INDEL_exonic"}{"unknown"}{$pid} = 0;
    $output_hash{$gene}{"INDEL_intronic"}{"."}{$pid} = 0;
    $output_hash{$gene}{"INDEL_ncRNA_exonic"}{"."}{$pid} = 0;
    $output_hash{$gene}{"INDEL_ncRNA_intronic"}{"."}{$pid} = 0;
    $output_hash{$gene}{"INDEL_ncRNA_UTR5"}{"."}{$pid} = 0;
    $output_hash{$gene}{"INDEL_splicing"}{"."}{$pid} = 0;
    $output_hash{$gene}{"INDEL_upstream"}{"."}{$pid} = 0;
    $output_hash{$gene}{"INDEL_downstream"}{"."}{$pid} = 0;
    $output_hash{$gene}{"INDEL_UTR3"}{"."}{$pid} = 0;
    $output_hash{$gene}{"INDEL_UTR5"}{"."}{$pid} = 0;
}

sub setup_gene_indel_germline {
    my $gene = shift;
    my $pid = shift;
    $output_hash_germline{$gene}{"INDEL_exonic"}{"frameshift deletion"}{$pid} = 0;
    $output_hash_germline{$gene}{"INDEL_exonic"}{"frameshift insertion"}{$pid} = 0;
    $output_hash_germline{$gene}{"INDEL_exonic"}{"nonframeshift deletion"}{$pid} = 0;
    $output_hash_germline{$gene}{"INDEL_exonic"}{"nonframeshift insertion"}{$pid} = 0;
    $output_hash_germline{$gene}{"INDEL_exonic"}{"stopgain"}{$pid} = 0;
    $output_hash_germline{$gene}{"INDEL_exonic"}{"unknown"}{$pid} = 0;
    $output_hash_germline{$gene}{"INDEL_intronic"}{"."}{$pid} = 0;
    $output_hash_germline{$gene}{"INDEL_ncRNA_exonic"}{"."}{$pid} = 0;
    $output_hash_germline{$gene}{"INDEL_ncRNA_intronic"}{"."}{$pid} = 0;
    $output_hash_germline{$gene}{"INDEL_ncRNA_UTR5"}{"."}{$pid} = 0;
    $output_hash_germline{$gene}{"INDEL_splicing"}{"."}{$pid} = 0;
    $output_hash_germline{$gene}{"INDEL_upstream"}{"."}{$pid} = 0;
    $output_hash_germline{$gene}{"INDEL_downstream"}{"."}{$pid} = 0;
    $output_hash_germline{$gene}{"INDEL_UTR3"}{"."}{$pid} = 0;
    $output_hash_germline{$gene}{"INDEL_UTR5"}{"."}{$pid} = 0;
}
