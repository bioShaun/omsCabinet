#!/usr/bin/perl

# add strand information to repeat bed file 
# 2018.5.15

use Bio::EnsEMBL::Registry;
use Getopt::Long;

GetOptions(
    'species|s=s' => \$species,
    'out|o=s'     => \$out_dir,
    'help|h!'     => \$help,
    );

if($help){
    print "###################################################################\n";
    print "usage perl $0 -s species -o out_dir\n";
    exit(1)
}
    
my $registry = 'Bio::EnsEMBL::Registry';
my %strand_hash = (
    1 => '+',
    0 => '.',
    -1 => '-',
    );

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);

$registry->set_reconnect_when_lost(1);

my $output="${out_dir}/${species}.repeat.bed3";
my $process="${out_dir}/downloaded.items";

if ( ! -e ${out_dir}) {
    mkdir( $out_dir ) or die "无法创建 $out_dir 目录, $!";
}

my %exists_item;
if(-e $process) {
    open(PROC, "<$process");
    while(<PROC>){
	chomp;
	$exists_item{$_} = 1;
    }
    close(PROC) || die "无法关闭文件";
}

if ( exists $exists_item{"FINISHED"} ) {
    printf ("%s:download finished.\n", $species);
}
else {
    my $slice_adaptor = $registry->get_adaptor( $species, 'Core', 'Slice' );
    my $all_chrom = $slice_adaptor->fetch_all('chromosome');
    my $all_scaffold = $slice_adaptor->fetch_all('scaffold');
    
    my @all_chrom = @{$all_chrom};
    my @all_scaffold = @{$all_scaffold};
    my @all_items = (@all_chrom, @all_scaffold);
    
    open(DATA, "+>>$output");
    open(PROC, "+>>$process");
    foreach my $chrom (@all_items) {
        # Obtain a slice covering the entire chromosome X
	my $chrom_id = $chrom->seq_region_name();
        next if exists $exists_item{"$chrom_id"};
        my @repeats = @{ $chrom->get_all_RepeatFeatures() };
        foreach my $repeat (@repeats) {
            my $strand = $repeat->strand();
            my $out_strand = $strand_hash{$strand};
            my $consensus = $repeat->repeat_consensus();
            my $repeat_name = $consensus->name();
            my $repeat_class = $consensus->repeat_class();
            my $repeat_type = $consensus->repeat_type(); 
            printf DATA ( "%s\t%d\t%d\t%s\t%s\t%s\t%s\n",
    		$chrom_id, $repeat->start(), $repeat->end(), 
                $out_strand ,$repeat_name, $repeat_class, $repeat_type);
        }
        printf PROC "$chrom_id\n";
    }
    
    printf PROC "FINISHED\n";
    close(DATA) || die "无法关闭文件";
    close(PROC) || die "无法关闭文件";
    
}


