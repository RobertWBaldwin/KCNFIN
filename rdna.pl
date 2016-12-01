#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use Graph; 
use RDseq;
use Carp; 
use Getopt::Long;
use File::Basename;
use lib dirname(__FILE__).'/..';


use constant LINK => "n" x 50;
no warnings "experimental::autoderef";
my $REPORT;
my $REPORT_FH;
my $tmpDIR;
my $outDIR;
my %ALIGN;
my %RNA_SEQ;
my %DNA_SEQ; 
my $K; 
my $C;

my ( $DNA_FILE, $RNA_FILE, $PSL_FILE, $PATH_SCORE_MIN, $NODE_SCORE_MIN, $PID_MIN, $DNA_MAX, $GAP_PID_MAX, $GAP_DNA_MAX, $INTRON_MAX, $OVERLAP, $scaffold, $gap, $redundancy, $help ); 
# default values for scaffold options
my $DEFAULT_TMP_DIR;
my $DEFAULT_OUT_DIR; 
my $DEFAULT_PATH_SCORE = 50.1;
my $DEFAULT_NODE_SCORE = 1.6;
my $DEFAULT_PID        = 98.1;
my $DEFAULT_DNA_MAX    = 50;
my $DEFAULT_RNA_MAX    = 100; 
my $DEFAULT_INTRON_MAX = 100000;
# default values for gap options
my $DEFAULT_GAP_PID    = 98.0; 
my $DEFAULT_GAP_DNA    = 15;
my $DEFAULT_OVERLAP    = 5;
my $match =  101;

# remember idea for having different scores ( coverage score, low penalty alignment score, high penalty alignment score ) allowing users to select what scoring system works best for their data ( i.e., no one scoring system may work the best for all input ) 



GetOptions
(
    "dna=s"      => \$DNA_FILE,
    "rna=s"      => \$RNA_FILE,
    "psl=s"      => \$PSL_FILE,
    "scaffold"   => \$scaffold,
    "gap"        => \$gap,
    "redundancy" => \$redundancy,
    "pscore=f"   => \$PATH_SCORE_MIN,
    "nscore=f"   => \$NODE_SCORE_MIN,
    "spid=f"     => \$PID_MIN,
    "sdcount=i"  => \$DNA_MAX,
    "apid=f"     => \$GAP_PID_MAX,
    "adcount=i"  => \$GAP_DNA_MAX,
    "intron=i"   => \$INTRON_MAX,
    "overlap=i"  => \$OVERLAP,
    "h|help"     => \$help,
);


print <<'MESSAGE';

RDNA_Scaffolder: A tool for scaffolding genome assemblies using RNA sequence data

Author: Robert Baldwin ( rb14sp@brocku.ca )

MESSAGE


my $USAGE = <<'USAGE';

    Usage:  --scaffold  --psl <string>  --rna <string>  --dna <string>  [options] 
    
            --scaffold        Takes both DNA and RNA fasta file and regular psl file describing the alignments between the sequences in the given fasta files  
            --gap             Takes same input as -scaffold option and will attempt to fill gaps and update alignment information prior to scaffolding          
            --dna             DNA fasta file
            --rna             RNA fasta file 
            --psl             Standard BLAT psl desribing alignments between sequences located in DNA and RNA fasta files
            -h|--help         prints this message

  [options]:

           [--pscore <Num> scaffolding path coverage score used to filter paths, ranges from 0 to 100 ( default 50 )]
           [--nscore <Num> alignment coverage score used to filter alignments from psl, ranges from 0 to 100 ( default 1.6 )]   
           [--pid <Num>   percentage of identity used to filter psl alignments ( default 98.1 )]
           [--dcount <Int> maximum number of DNA sequences that each RNA sequence can be aligned to ( default 50 )]  
           [--intron <Int> maximum intron size that can be represented by an edge ( default 100000 )]  
           [--overlap <Int> maximum base pairs that two DNA sequences can be overlapping on RNA before overlap is considered significant ( default 5 )] 

    NOTE: requires blat and bl2seq be in default path 

USAGE


die $USAGE if ! $scaffold and ! $gap and ! $help;
 
if ( ! $DNA_FILE or ! $RNA_FILE or ! $PSL_FILE ) 
{
    print qq[FILES!\n];
    die $USAGE;
}

croak qq[Can't find DNA fasta file: $DNA_FILE\n] unless -f $DNA_FILE || -l $DNA_FILE; 
croak qq[Can't find RNA fasta file $RNA_FILE\n] unless -f $RNA_FILE || -l $RNA_FILE;
croak qq[Can't find PSL file $PSL_FILE\n] unless -f $PSL_FILE || -l $PSL_FILE;
$PATH_SCORE_MIN = defined $PATH_SCORE_MIN ? $PATH_SCORE_MIN : $DEFAULT_PATH_SCORE;
$NODE_SCORE_MIN = defined $NODE_SCORE_MIN ? $NODE_SCORE_MIN : $DEFAULT_NODE_SCORE;
$DNA_MAX = defined $DNA_MAX ? $DNA_MAX: $DEFAULT_DNA_MAX;
$PID_MIN = defined $PID_MIN ? $PID_MIN : $DEFAULT_PID;  
$INTRON_MAX = defined $INTRON_MAX ? $INTRON_MAX : $DEFAULT_INTRON_MAX;
$OVERLAP = defined $OVERLAP ? $OVERLAP: $DEFAULT_OVERLAP;

if ( defined $gap )
{
    $GAP_PID_MAX = defined $GAP_PID_MAX ? $GAP_PID_MAX : $DEFAULT_GAP_PID;
    $GAP_DNA_MAX = defined $GAP_DNA_MAX ? $GAP_DNA_MAX : $DEFAULT_GAP_DNA;
}

my $blat = qw[/work/lianglab/bin/blat];
open $REPORT_FH, '>', \$REPORT;
$tmpDIR = defined $tmpDIR ? $tmpDIR/"tmp_$$" : "tmp_$$";
$outDIR = defined $outDIR ? $outDIR/"RDNA_$$" : "RDNA_$$"; 
mkdir $tmpDIR, or die "cannot create $tmpDIR: $!" unless -e $tmpDIR;
mkdir $outDIR, or die "cannot create $outDIR: $!" unless -e $outDIR;
print qq[Writing temporary RNA file.\n]; 
inFasta ( $RNA_FILE, \%RNA_SEQ );
open my ($rnafh), '>', qq[$tmpDIR/RNA_tmp.fa] or die "$!\n";
while( my ( $id, $seq ) = each %RNA_SEQ )
{
   # $seq = formatFasta( \$seq ); 
    print $rnafh ">$id\n$seq\n";
}
close $rnafh;
print qq[Writing temporary DNA file.\n];
inFasta ( $DNA_FILE, \%DNA_SEQ );
open my ($dnafh), '>', qq[$tmpDIR/DNA_tmp.fa] or die "$!\n";
while ( my ( $id, $seq ) = each %DNA_SEQ )
{
    print $dnafh ">$id\n$seq\n";
}
close $dnafh;
$C =()= keys %RNA_SEQ;
print qq[$C rna in $RNA_FILE.\n];
$C = () = keys %DNA_SEQ;
print qq[$C dna in $DNA_FILE.\n];  
print qq[Retreiving alignment information.\n];
inPSL ( \%ALIGN, $PSL_FILE, $NODE_SCORE_MIN, $PID_MIN, $DNA_MAX );
print qq[$C DNA from psl\n]; 
$RNA_SEQ{$_} = () foreach keys %RNA_SEQ;


=some
for my $dna ( keys %{$ALIGN{$rna}} )
{
    print qq[$dna\n];
    print qq[Matches: $ALIGN{$rna}{$dna}->{matches}\nMismatches: $ALIGN{$rna}{$dna}->{mismatches}\nRepmatches: $ALIGN{$rna}{$dna}->{repmatches}\nNcount: $ALIGN{$rna}{$dna}->{ncount}\nQinsert: $ALIGN{$rna}{$dna}->{qinsert}\tQbase: $ALIGN{$rna}{$dna}->{qbase}\nTinsert: $ALIGN{$rna}{$dna}->{tinsert}\tTbase: $ALIGN{$rna}{$dna}->{tbase}\nStrand: $ALIGN{$rna}{$dna}->{strand}\nBlocknum: $ALIGN{$rna}{$dna}->{bnum}\nQstart: $ALIGN{$rna}{$dna}->{qstart}\tQend: $ALIGN{$rna}{$dna}->{qend}\nTstart: $ALIGN{$rna}{$dna}->{tstart}\tTend: $ALIGN{$rna}{$dna}->{tend}\nRNA length: $ALIGN{$rna}{$dna}->{qsize}\nDNA length: $ALIGN{$rna}{$dna}->{tsize}\nPercent Coverage: $ALIGN{$rna}{$dna}->{pcov}\n];
    print qq[$_\t] foreach @{$ALIGN{$rna}{$dna}->{q_starts}};
    print qq[\n];
    print qq[$_\t] foreach @{$ALIGN{$rna}{$dna}->{t_starts}};
    print qq[\n];
    print qq[$_\t] foreach @{$ALIGN{$rna}{$dna}->{b_sizes}}; 
    print qq[\n\n];


    foreach ( @{$ALIGN{$rna}{$dna}->{blocks}} )
    {
        print qq[RNA Start: $_->{qstart}\tRNA end: $_->{qend}\nDNA Start: $_->{tstart}\tDNA End: $_->{tend}\n\n];
        my $rna_exon = substr $RNA_SEQ{$rna}, $_->{qstart} - 1, $_->{qend} - ( $_->{qstart} - 1 );
        print qq[RNA exon:\n$rna_exon\n];
        my $seq = $ALIGN{$rna}{$dna}->{strand} eq '+' ? $DNA_SEQ{$dna} : rc ( $DNA_SEQ{$dna} );
        my $dna_exon = substr $seq, $_->{tstart} - 1, $_->{tend} - ( $_->{tstart} - 1 ); 
        print qq[DNA coding seq:\n$dna_exon\n];
    }
    print qq[\n\n];
}
=cut



if ( $C == 0 )
{
    print qq[0 alignments found.\n];
    exit;
}

for my $rna ( keys %ALIGN )
{
    $K->{$_}++ foreach ( keys %{$ALIGN{$rna}} );
}
$C =()= keys %{$K};
print qq[$C dna in alignments\n];
print qq[removing dna sequences with zero alignments\n]; 
for my $dna ( keys %DNA_SEQ )
{
    delete $DNA_SEQ{$dna} if ! exists $K->{$dna};
}
$C =()= keys %DNA_SEQ;
print qq[$C dna sequences remaining\n];
$K = ();
print qq[removing rna sequences with zero alignments\n];
for my $rna ( keys %RNA_SEQ )
{
    delete $RNA_SEQ{$rna} if ! exists $ALIGN{$rna};
}
$C =()= keys %RNA_SEQ;
print qq[$C rna sequences remaining.\n];
print qq[Grouping rna/dna sequence elements based on alignment data\n];
my $groups_ref = groupSequences ( \%ALIGN );
$C = scalar @$groups_ref;
print qq[$C groups created\n];
print qq[Filtering groups\n];
my @groups;

for ( my $i = 0; $i < scalar @$groups_ref; $i++ )
{
    my $group = @$groups_ref[$i];
    my @RNA = @{$group->{rna}};
    my @DNA = @{$group->{dna}};

    if ( scalar @DNA == 1  ) 
    {
        delete $DNA_SEQ{$_} foreach @DNA;
        
        foreach ( @RNA )
        {
            delete $ALIGN{$_};
            delete $RNA_SEQ{$_};
        } 
        next; 
    }

    push @groups, $group;  
}  

$C = @groups;
print qq[$C groups remaining\n];
$C =()= keys %DNA_SEQ;
print qq[$C dna remaining\n];
$C =()= keys %RNA_SEQ;
print qq[$C rna remaining\n];
print qq[creating objects for remaining sequence data\n]; 

for my $rna ( keys %ALIGN ) 
{
    $RNA_SEQ{$rna} = RNA->new ( status => 'r' );

    for my $dna ( keys %{$ALIGN{$rna}} ) 
    {
        push @{$RNA_SEQ{$rna}->{alignments}}, $ALIGN{$rna}{$dna};
        $ALIGN{$rna}{$dna} = 1;
        push @{$K->{$dna}}, $rna; 
    }
}


for my $dna ( keys %DNA_SEQ )
{
    ( my $rc = reverse $DNA_SEQ{$dna} ) =~ tr /ATCGatcg/TAGCtagc/;
    $DNA_SEQ{$dna} = DNA->new ( seq => $DNA_SEQ{$dna}, revcom => $rc,  length => length $DNA_SEQ{$dna}, rna => \@{$K->{$dna}}, gap => $DNA_SEQ{$dna} =~ tr/N//, status => 'o' );  
    push @{$DNA_SEQ{$dna}->{predecessors}}, $dna;   
} 


=some
for my $rna ( keys %RNA_SEQ )
{
    for my $alignment ( @{$RNA_SEQ{$rna}->{alignments}} )
    {
        if ( $alignment->{low_score} < 3 )
        {
            my $dna = $alignment->{dna};

            if ( scalar @{$DNA_SEQ{$dna}->{rna}} > 1 || $alignment->{pid} < 99.0  )
            {
                @{$DNA_SEQ{$dna}->{rna}} = grep { $_ ne $rna } @{$DNA_SEQ{$dna}->{rna}};
                @{$RNA_SEQ{$rna}->{alignments}} = grep { $_->{dna} ne $dna } @{$RNA_SEQ{$rna}->{alignments}};   
                delete $ALIGN{$rna}{$dna};
            }          
        }
    }
}
=cut


####### gap filling and assembly part starts here. problem is too slow for large input because all against all alignment with bl2seq takes long time 
####### would be faster to reduce number of alignments by only trying sequences that overlap on the RNA 
####### you can't update alignment information as you go along so any updated sequence would need to inherit alignment lists of it predecessors
####### worth trying....

#my $start_time = time;

=some

for ( my $i = 0; $i < scalar @groups; $i++ )
{
    my $group = $groups[$i];
    my @RNA = @{$group->{rna}};
    my $id = idGenerate ( 'RDNAU', $i );
    my $seen; 

    for my $rna ( @RNA )
    {
        my @DNA = keys %{$ALIGN{$rna}};
        $C = @DNA;
        next unless $C > 1; 
        print qq[Processing $C sequences for $rna in group $i\n];
        #my $x_graph = Graph::Undirected->new();
        #my %A = map { $_->{dna} => { start => $_->{qstart}, end => $_->{qend} } } @{$RNA_SEQ{$rna}->{alignments}}; 

        for my $dna ( @DNA )
        {
            next if $DNA_SEQ{$dna}->{file};  
            open my ( $os_fh ), '>', qq[$tmpDIR/$dna.fa] or die "$!\n";
            print $os_fh qq[>$dna\n$DNA_SEQ{$dna}->{seq}\n];  
            close $os_fh; 
            $DNA_SEQ{$dna}->{file} = 1;
            #my @overlap = grep { $A{$_} ne $dna && ! ( $A{$_}->{start} > $A{$dna}->{end} || $A{$_}->{end} < $A{$dna}->{start} ) } keys %      
        }

        RDassemble ( \@DNA, $seen, $id );
    } 
}    

$C = 0;
$K = 0;

for my $dna ( keys %DNA_SEQ )
{
    $C++ if $DNA_SEQ{$dna}->{status} eq 'ur' || $DNA_SEQ{$dna}->{status} eq 'uu' || $DNA_SEQ{$dna}->{status} eq 'u';
    $K++ if $DNA_SEQ{$dna}->{status} eq 'o'; 
}

print qq[$K original sequences remaining\n];
print qq[$C updated sequences created\n]; 
$C = 0;
$K = 0;

for my $dna ( keys %DNA_SEQ )
{
    $C++ if $DNA_SEQ{$dna}->{status} eq 'u';
}
# updated sequences not subsequently removed either because they were redundant or assembled into larger sequence
print qq[$C primary updated sequences\n]; 
$C = 0; 

for my $dna ( keys %DNA_SEQ )
{
    $C++ if $DNA_SEQ{$dna}->{status} eq 'or';
    $K++ if $DNA_SEQ{$dna}->{status} eq 'ur';
}

print qq[$C original and $K updated sequences redundant\n];
$C = 0;

# Determine how many alignments need to be updated 
for my $rna ( keys %ALIGN )
{
    for my $dna ( keys %{$ALIGN{$rna}} )
    {
        $C++ if $ALIGN{$rna}{$dna} eq 'update'; 
    }
}

$K = ();
# Obtain the RNA sequence data, keeping only those sequences in alignment data 
$RNA_FILE = qq[$tmpDIR/RNA_tmp.fa]; 
inFasta ( $RNA_FILE, \%{$K} );

for my $rna ( keys %{$K} )
{
    if ( ! exists $ALIGN{$rna} )
    {
        delete $K->{$rna};
    }

    else 
    {
        $RNA_SEQ{$rna}->{seq} = $K->{$rna}; 
    }     
}

$K = ();

# If we have updated sequences ( and so alignments to update )
if ( $C )
{ 
    my $update_dna = qq[$tmpDIR/updated_sequences.fa];
    print qq[$C alignments to be updated\nprinting updated dna sequences to $update_dna\n];
    open my ($update_fh), '>', qq[$update_dna] or die "$!\n";
    $C = 0; 
    # print updated sequences to file
    for my $dna ( keys %DNA_SEQ )
    {
        if ( $DNA_SEQ{$dna}->{status} eq 'u' )
        {
            print $update_fh qq[>$dna\n$DNA_SEQ{$dna}->{seq}\n];
            $K->{$_} = () foreach @{$DNA_SEQ{$dna}->{rna}};
        }
          
    }

    # print RNA to file
    $RNA_FILE = qq[$tmpDIR/RNA_new.fa];
    open my ($rna_fh), '>', qq[$RNA_FILE] or die "$!\n";
    print $rna_fh qq[>$_\n$RNA_SEQ{$_}->{seq}\n] foreach keys %{$K};   
    close $rna_fh;
    $K = ();
    print qq[Updating alignment information\nCalling BLAT with $RNA_FILE and $DNA_FILE\n];
    my $psl_out = 'update.psl'; 
    system ( $blat, '-t=dna', '-q=rna', '-fine', "$DNA_FILE", "$RNA_FILE", "$psl_out" );
    $C = inPSL ( \%{$K}, $psl_out, $NODE_SCORE_MIN, $PID_MIN, $DNA_MAX );
    print qq[NO ALIGNMENTS\n] if ! $C;
    $C = 0; 

    for my $rna ( keys %ALIGN )
    {
        for my $dna ( keys %{$ALIGN{$rna}} )
        {
            if ( $ALIGN{$rna}{$dna} eq 'update' )
            {
                if ( exists $K->{$rna}->{$dna} )
                {
                    push @{$RNA_SEQ{$rna}->{alignments}}, $K->{$rna}->{$dna}; 
                    $ALIGN{$rna}{$dna} = 1; 
                }

                else
                {
                    @{$DNA_SEQ{$dna}->{rna}} = grep { $_ ne $rna } @{$DNA_SEQ{$dna}->{rna}};
                    delete $ALIGN{$rna}{$dna};
                    $C++; 
                }
            }
        }
    }
}

$K = ();

if ( $C )
{
    print qq[$C alignments not updated\n];
    open my ($dna_fh), '>', qq[$outDIR/excluded_sequences.fa] or die "$!\n";
    
    # look for updated dna sequences with no alignments, print them to file and restore all original sequences         
    for my $dna ( keys %DNA_SEQ )
    {   
        my @rna = @{$DNA_SEQ{$dna}->{rna}}; 
        my $c = scalar @rna;
        # If updated sequence does not have any alignment with RNA ( alignment quality went below thresehold )
        if ( ! $c )
        {
            print qq[$dna being restored to original sequences, no RNA alignment\n];
 
            if ( $DNA_SEQ{$dna}->{status} eq 'u' )
            { 
                print $dna_fh qq[>$dna | @{$DNA_SEQ{$dna}->{predecessors}}\n$DNA_SEQ{$dna}->{seq}\n];
                $DNA_SEQ{$dna}->{status} = 'x';
                # restore all the original sequences used to create the updated sequence
                foreach ( @{$DNA_SEQ{$dna}->{predecessors}} )
                {
                    $DNA_SEQ{$_}->{status} = 'o'; 
                    $ALIGN{$_}{$dna} = 1 foreach @rna;  
                }
            }

            elsif ( $DNA_SEQ{$dna}->{status} eq 'o' )
            {
                print qq[Error: original sequence with no alignments\n]; 
            }
        }
    }
}

my $end_time = time;
my $total_time = $end_time - $start_time;
print qq[TIME: $total_time\n];
$C =()= keys %ALIGN;   
print qq[$C rna remaining\n]; 
print qq[Regrouping sequences\n];
$groups_ref = groupSequences ( \%ALIGN );
$C = @$groups_ref; 
print qq[$C groups created\n];
=cut




$K = ();
$C = ();
@groups = sort { @{$a->{rna}} <=> @{$b->{rna}} } @groups; 

for ( my $i = 0; $i < scalar @groups; $i++ )
{
    my $group = $groups[$i];
    my @RNA = @{$group->{rna}};
    my $scaffold_id = idGenerate('KCNFIN', $i); 
    my @E = ();
    $K = scalar @RNA; 
    $C = (); 
    $K = (); 
    # determine all possible scaffolding edges for each RNA      
    for my $rna ( @RNA )
    {   
        my @dna = keys %{$ALIGN{$rna}};
        next if scalar @dna == 1;
        @{$RNA_SEQ{$rna}->{alignments}} = sort { $a->{qstart} <=> $b->{qstart} } @{$RNA_SEQ{$rna}->{alignments}};
      
        for ( my $i = 0; $i < scalar @{$RNA_SEQ{$rna}->{alignments}} - 1; $i++ )
        {
            my $alignment_a = @{$RNA_SEQ{$rna}->{alignments}}[$i];

            for ( my $j = $i + 1; $j < scalar @{$RNA_SEQ{$rna}->{alignments}}; $j++ )
            {
                my $alignment_b = @{$RNA_SEQ{$rna}->{alignments}}[$j]; 
                # don't add edge if there's significant overlap on RNA
                next if $alignment_a->{qend} > $alignment_b->{qstart} + $OVERLAP;
                # don't add edge if there's a third DNA sequence aligned to the RNA between the two 
                next if grep { $_->{qstart} >= $alignment_a->{qend} && $_->{qend} <= $alignment_b->{qstart} } @{$RNA_SEQ{$rna}->{alignments}}; 
                push @E, { v1 => $alignment_a->{dna}, v2 => $alignment_b->{dna}, rna => $rna, state => 1 };            
            }
        }
    }       

    @RNA = ();   
    next if ! scalar @E;

    for my $edge ( @E )
    {
        my $tend;
        my $tstart;
        my $qstart;
        my $qend;
       
        for my $a ( @{$RNA_SEQ{$edge->{rna}}->{alignments}} )
        {
            if ( $a->{dna} eq $edge->{v1} )
            { 
		$tend = $a->{tend};
                $qend = $a->{qend};
            }

            if ( $a->{dna} eq $edge->{v2} )
            {
                $tstart = $a->{tstart};
                $qstart = $a->{qstart}; 
            }
            
        }
 
        # remove large introns only for edges where edge doesn't span a gap in the alignment with the RNA ( edge needs to actually represent an intron )  
        next unless $qstart <= $qend;           
        my $min_intron = ( $DNA_SEQ{$edge->{v1}}->{length} - $tend ) + $tstart;
        # print qq[Removed edge $edge->{v1} --> $edge->{v2}\t$min_intron\n] if $min_intron > $INTRON_MAX; 
        $edge->{state} = 0 if $min_intron > $INTRON_MAX;   
    }
 
    @E = grep { $_->{state} } @E;
    @RNA = keys { map { $_->{rna} => 1 } @E }; 
    next if ! scalar @E;        
    my @MAX_PATHS;     

    # For each RNA select the maximum scoring path: graph all possible edges; determine all possible paths; select highest scoring path  
    for my $rna ( @RNA )
    {
        my $g = Graph->new(); 
        my @e = grep { $_->{rna} eq $rna } @E;
        my %weights = map { $_->{dna} => $_->{low_score} } @{$RNA_SEQ{$rna}->{alignments}};
        my %matches = map { $_->{dna} => $_->{matches} } @{$RNA_SEQ{$rna}->{alignments}};

        for my $edge ( @e )
        {
            $g->add_weighted_vertices( $edge->{v1}, $weights{$edge->{v1}}, $edge->{v2}, $weights{$edge->{v2}} );
            $g->add_edge ( $edge->{v1}, $edge->{v2} );
        }
        
        my @source = $g->source_vertices;
        my @paths; 
        # determine all possible scaffolding paths       
        for my $source_node ( @source )
        {
           my $dfs_paths = DFS( $source_node, $g );
           push @paths, $_ foreach @$dfs_paths; 
        }

        # remove any paths whose set of vertices is present in a longer path 
        my $filt_paths = ( scalar @paths > 1 ) ? set_inclusion(\@paths) : \@paths;  
        my @path_score; 
        
        # score the paths
        for my $path ( @$filt_paths )
        {
            my $score;
            $score += $g->get_vertex_weight($_) foreach @$path;  
            my $m;
            $m += $matches{$_} foreach @$path;
            push @path_score, { rna => $rna, score => $score, path => $path, matches => $m  };                 
        } 

        # select the top scoring path 
        @path_score = sort { $b->{score} <=> $a->{score} } @path_score;
        my $max_path_a = $path_score[0];
        my $max_path_b = exists $path_score[1] ? $path_score[1] : undef;
        @{$RNA_SEQ{$rna}->{alignments}} = sort { $b->{low_score} <=> $a->{low_score} } @{$RNA_SEQ{$rna}->{alignments}};
        # skip RNA whose max scoring path is too low
        if ( $max_path_a->{score} < $PATH_SCORE_MIN )
        {
            $RNA_SEQ{$rna}->{status} = 's';
            next;
        }
        # skip RNA whose max scoring path is less than the top scoring alignment + minimum node score     
        next if $max_path_a->{score} < ${$RNA_SEQ{$rna}->{alignments}}[0]->{low_score} + $NODE_SCORE_MIN;
      
=some 
        if ( defined $max_path_b )
        {
            my $score_difference = $max_path_a->{score} - $max_path_b->{score};
            my $match_difference = $max_path_a->{matches} - $max_path_b->{matches};

            if ( $score_difference < 2 || $match_difference < 100 )
            {
                my %scores = map { $_->{dna} => $_->{coverage_score} } @{$RNA_SEQ{$max_path_a->{rna}}->{alignments}};
                my %pid = map { $_->{dna} => $_->{pid} } @{$RNA_SEQ{$max_path_a->{rna}}->{alignments}}; 
                my %q = map { $_->{dna} => $_->{qinsert} } @{$RNA_SEQ{$max_path_a->{rna}}->{alignments}};  
                my %m = map { $_->{dna} => $_->{matches} } @{$RNA_SEQ{$max_path_a->{rna}}->{alignments}};   
               # print qq[$rna\nMax A: $max_path_a->{score}\n];
               # print qq[$_ $scores{$_} $pid{$_} $q{$_} $m{$_}\t] foreach @{$max_path_a->{path}};
               # print qq[\nMax B: $max_path_b->{score}\n];
               # print qq[$_ $scores{$_} $pid{$_} $q{$_} $m{$_}\t] foreach @{$max_path_b->{path}};  
               # print qq[\n\n];
            }
        }
=cut

        my $path = PATH->new ( rna => $rna, score => $max_path_a->{score}  );
        my %strands = map { $_->{dna} => { str => $_->{strand}, scr => $_->{low_score} } } @{$RNA_SEQ{$rna}->{alignments}};
        push @{$path->{path}}, { id => $_, strand => $strands{$_}->{str}, score => $strands{$_}->{scr} } foreach @{$max_path_a->{path}};
        # collect all the max paths for each group 
        push @MAX_PATHS, $path; 
    }

    @E = ();
    @RNA = ();
    next if ! scalar @MAX_PATHS;

    # group has one max path   
    if ( scalar @MAX_PATHS == 1 )
    {
        scaffold_a( $scaffold_id, $MAX_PATHS[0] );
        next;
    }
 
    # group has multiple max paths; do some share vertices? 
    my $undir = Graph::Undirected->new();

    for my $path ( @MAX_PATHS )
    {
        for ( my $i = 0; $i < scalar @{$path->{path}} - 1; $i++ )
        {
            $undir->add_edge ( @{$path->{path}}[$i]->{id}, @{$path->{path}}[$i+1]->{id} );  
	}
    }	 
             
    my @cc = $undir->connected_components();
    my @edges = $undir->edges; 
    $undir->delete_edge ( @$_[0], @$_[1] ) foreach @edges; 
        
    for my $component ( @cc )
    {   
        my %vertices;
        $vertices{$_} = () foreach @$component;
        my @component_paths;
        $K = ();

        for my $path ( @MAX_PATHS )
        {            
            for my $dna ( @{$path->{path}} )
            {
                if ( exists $vertices{$dna->{id}} )
                {
                    push @component_paths, $path;
                    last;
                }     
            }
        }         

        # max paths may not share vertices; some components have a single max path 
        if ( scalar @component_paths == 1 )
        {
            scaffold_a( $scaffold_id, $component_paths[0] );
            next;
        }
                           
        # components here are "complex" consisting of multiple max paths 
        # first step is to find the most plausible order for these paths based on the strand positions of the dna sequences that they have in common (+/-)
        my %plus;
        my %minus;  
        plus_minus ( \@component_paths, \%plus, \%minus );           
        my @plus_minus = grep { exists $plus{$_} && exists $minus{$_} } keys %vertices;          

        if ( scalar @plus_minus )
        {
            print qq[Possible reverse/foreward transcripts\n];
            order_paths ( \@component_paths );
            %plus = ();
            %minus = ();            
            plus_minus ( \@component_paths, \%plus, \%minus ); 
            @plus_minus = grep { exists $plus{$_} && exists $minus{$_} } keys %vertices;
 
            if ( scalar @plus_minus )  
	    {
                print qq[Plus Minus not resolved\n];
                next; 
	    }   
        }

        # graph components that have a plausible path order  
        my $digraph = Graph->new(); 
        my @foreward_paths = grep { $_->{order} eq 'foreward' } @component_paths;
        my @reverse_paths = grep { $_->{order} eq 'reverse' } @component_paths;

        for my $path ( @foreward_paths )
        {
            for ( my $i = 0; $i < scalar @{$path->{path}} - 1; $i++ )
            {
                if ( ! $digraph->has_edge( @{$path->{path}}[$i]->{id}, @{$path->{path}}[$i+1]->{id} ) )
		{
                    $digraph->add_edge ( @{$path->{path}}[$i]->{id}, @{$path->{path}}[$i+1]->{id} );
		}
            }
        } 
	
        for my $path ( @reverse_paths )
        {
            my @reverse = reverse @{$path->{path}};

            for ( my $i = 0; $i < scalar @reverse - 1; $i++ )
            {
                if ( ! $digraph->has_edge( $reverse[$i]->{id}, $reverse[$i+1]->{id} ) )
		{
                    $digraph->add_edge ( $reverse[$i]->{id}, $reverse[$i+1]->{id} );
		}
            } 
        }

        if ( $digraph->is_cyclic )
        {
           print qq[Graph is Cyclic ################################\n];
           next; 
        }

        # get consensus path, if it exists
        my @raw_paths;    
        my @source = $digraph->source_vertices();
                
        for my $source_vertex ( @source )
        {
            my $dfs = DFS ( $source_vertex, $digraph );
            push @raw_paths, $_ foreach @$dfs; 
        }

        my $filt_paths = ( scalar @raw_paths > 1 ) ? set_inclusion( \@raw_paths ) : \@raw_paths;     
        @$filt_paths = sort { $a->[0] cmp $b->[0] || $a->[1] cmp $b->[1] || $a->[2] cmp $b->[2] || $a->[3] cmp $b->[3] } @$filt_paths;
  
        if ( scalar @$filt_paths == 1 )
        {                  
            scaffold_b( $scaffold_id, $filt_paths, \%plus, \%minus );
            next;
        }

        ### No consensus path; fall back on heuristics 

        $K = ();
        @raw_paths = ();
        @source = (); 

        # get the edge counts 
        for my $path ( @component_paths )
        {
            for ( my $i = 0; $i < scalar @{$path->{path}} - 1; $i++ )
            {
                if ( $path->{order} eq 'foreward' )
                {
                    $K->{@{$path->{path}}[$i]->{id}}->{@{$path->{path}}[$i+1]->{id}}++;
                }
      
                else
                {
                    $K->{@{$path->{path}}[$i+1]->{id}}->{@{$path->{path}}[$i]->{id}}++;
                }
            }
        }
       
        for my $A ( keys %{$K} )
        {
            for my $B ( keys %{$K->{$A}} )
            {
                push @source, [$A, $B] if $K->{$A}->{$B} > 1;
            }
        }

        for my $path ( @component_paths )
        {
            my %pid = map {  $_->{dna} => $_->{pid} } @{$RNA_SEQ{$path->{rna}}->{alignments}};
            my %match = map { $_->{dna} => $_->{matches} } @{$RNA_SEQ{$path->{rna}}->{alignments}}; 
            print qq[\n$path->{rna}\t$path->{score}\t$path->{order}\n]; 
            print qq[$_->{id} $_->{strand} $_->{score} $pid{$_->{id}} $match{$_->{id}}\t] foreach @{$path->{path}}; 
            print qq[\n\n];
        }

        # graph the multi edges only and try again 
        if ( scalar @source )
        {
            my $final_graph = Graph->new();
            @source = ();
 
            my @foreward_paths = grep { $_->{order} eq 'foreward' } @component_paths;
            my @reverse_paths = grep { $_->{order} eq 'reverse' } @component_paths;

            for my $path ( @foreward_paths )
            {
                for ( my $i = 0; $i < scalar @{$path->{path}} - 1; $i++ )
                {
                    my $A = @{$path->{path}}[$i]->{id};
                    my $B = @{$path->{path}}[$i+1]->{id};

                    if ( ! $final_graph->has_edge($A, $B ) )
	  	    {
                        $final_graph->add_edge ( $A, $B ) if $K->{$A}->{$B} > 1;
		    }
                }
            } 
	
            for my $path ( @reverse_paths )
            {
                my @reverse = reverse @{$path->{path}};

                for ( my $i = 0; $i < scalar @reverse - 1; $i++ )
                {
                    my $A = $reverse[$i]->{id};
                    my $B = $reverse[$i+1]->{id};

                    if ( ! $final_graph->has_edge( $A, $B ) )
		    {
                        $final_graph->add_edge ( $A, $B ) if $K->{$A}->{$B} > 1; 
		    }
                } 
            }

         
            if ( $final_graph->is_cyclic )
            {
                print qq[Graph is Cyclic ################################\n];
                next; 
            }   
            # recompute the paths
            @source = $final_graph->source_vertices();

            for my $source_vertex ( @source )
            {
                my $dfs = DFS ( $source_vertex, $final_graph );
                push @raw_paths, $_ foreach @$dfs; 
            }
      
            my $final_paths = ( scalar @raw_paths > 1 ) ? set_inclusion( \@raw_paths ) : \@raw_paths;
            @$final_paths = sort { $a->[0] cmp $b->[0] || $a->[1] cmp $b->[1] || $a->[2] cmp $b->[2] || $a->[3] cmp $b->[3] } @$final_paths;

            for my $path ( @$final_paths )
            {
                print qq[$_\t] foreach @$path;
                print qq[\n];
            }    
      
            $K = ();
  
            for my $path ( @$final_paths )
            {
                $K->{$_}++ foreach @$path;
            }
            # if paths share dna select one path only otherwise scaffold all of them 
            if ( grep { $K->{$_} > 1 } keys %{$K} )
	    {
                @$final_paths = @$final_paths[0]; 
            }
           
            for my $scaffold_path ( @$final_paths ) 
            { 
                print qq[Scaffolding Path E\n];
                print qq[$_\t] foreach @$scaffold_path;
                print qq[\n\n];
            }

            scaffold_b( $scaffold_id, $final_paths, \%plus, \%minus );
            next;       
                                                             
        }

        # otherwise scaffold as many strong connections as possible  
        else
        {
            $K = ();
 
            my $final_graph = Graph->new();
            @source = (); 
            my @foreward_paths = grep { $_->{order} eq 'foreward' } @component_paths;
            my @reverse_paths = grep { $_->{order} eq 'reverse' } @component_paths;

            for my $path ( @foreward_paths )
            {
                my %match = map { $_->{dna} => $_->{matches} } @{$RNA_SEQ{$path->{rna}}->{alignments}}; 

                for ( my $i = 0; $i < scalar @{$path->{path}} - 1; $i++ )
                {
                    my $A = @{$path->{path}}[$i]->{id};
                    my $B = @{$path->{path}}[$i+1]->{id};

                    if ( ! $final_graph->has_edge( $A, $B ) )
	  	    {
                        $final_graph->add_edge ( $A, $B ) if $match{$A} > 500 && $match{$B} > 500; 
		    }
                }
            } 
	
            for my $path ( @reverse_paths )
            {
                my %match = map { $_->{dna} => $_->{matches} } @{$RNA_SEQ{$path->{rna}}->{alignments}}; 
                my @reverse = reverse @{$path->{path}};
        
                for ( my $i = 0; $i < scalar @reverse - 1; $i++ )
                {
                    my $A = $reverse[$i]->{id};
                    my $B = $reverse[$i+1]->{id};

                    if ( ! $final_graph->has_edge( $A, $B ) )
		    {
                         $final_graph->add_edge ( $A, $B ) if $match{$A} > 500 && $match{$B} > 500;   
		    }
                } 
            }
            
             
            # recompute the paths
            @source = $final_graph->source_vertices();
            next unless scalar @source; 
 
            if ( $final_graph->is_cyclic )
            {
                print qq[Graph is Cyclic ################################\n];
                next; 
            }

            for my $source_vertex ( @source )
            {
                my $dfs = DFS ( $source_vertex, $final_graph );
                push @raw_paths, $_ foreach @$dfs; 
            }
      
            print qq[RAW PATHS:\n];

            for my $path ( @raw_paths )
            {
                print qq[$_\t] foreach @$path;
                print qq[\n];
            }

            my $final_paths = ( scalar @raw_paths > 1 ) ? set_inclusion( \@raw_paths ) : \@raw_paths;
            @$final_paths = sort { $a->[0] cmp $b->[0] || $a->[1] cmp $b->[1] || $a->[2] cmp $b->[2] || $a->[3] cmp $b->[3] } @$final_paths;
            print qq[FINAL PATHS###\n];

            for my $path ( @$final_paths )
            {
                print qq[$_\t] foreach @$path;
                print qq[\n];
            }
 
            print qq[\n];

             $K = ();
  
            for my $path ( @$final_paths )
            {
                $K->{$_}++ foreach @$path;
            }

            if ( grep { $K->{$_} > 1 } keys %{$K} )
	    {
                @$final_paths = @$final_paths[0]; 
            }

            scaffold_b ( $scaffold_id, $final_paths, \%plus, \%minus );           
	    next;
	}

        print qq[UNRESOLVED\n]; 	    
    }
}




$K = ();
$C = ();
undef %ALIGN;
print qq[\n\nsequence_count\t\tdefined_bases\t\tN_bases\t\t\tn_bases\t\tTotal_Bases\n\n];
# open original dna fasta
inFasta ( $DNA_FILE, \%{$K} );
# record sequence data summary 
datFasta ( $DNA_FILE, \%{$K}, \%{$C} );
print qq[$DNA_FILE\n$C->{$DNA_FILE}->{seqcount}\t\t$C->{$DNA_FILE}->{defbases}\t\t$C->{$DNA_FILE}->{nbases}\t\t$C->{$DNA_FILE}->{lbases}\t\t$C->{$DNA_FILE}->{totbases}\n\n];
# print remaining sequences to file; those not in alignment data or not touched by program ( scaffolds or updated sequence, marked 'o' for original )
$DNA_FILE = qq[$outDIR/remainder.fa];
open my ($remainder_fh), '>', $DNA_FILE or die "$!\n";
for my $dna ( sort keys %{$K} )
{
    # sequence not in any alignment    
    if ( ! exists $DNA_SEQ{$dna} )
    {
        print $remainder_fh qq[>$dna\n$K->{$dna}\n];
    }     
    # or if it was and wasn't touched 
    elsif ( $DNA_SEQ{$dna}->{status} eq 'o' )
    { 
        print $remainder_fh qq[>$dna\n$DNA_SEQ{$dna}->{seq}\n]
    }
}
close $remainder_fh;
undef %{$K};
# load remainder dna file
inFasta ( $DNA_FILE, \%{$K} );
# record sequence data summary
datFasta ( $DNA_FILE, \%{$K}, \%{$C} );
print qq[$DNA_FILE\n$C->{$DNA_FILE}->{seqcount}\t\t$C->{$DNA_FILE}->{defbases}\t\t$C->{$DNA_FILE}->{nbases}\t\t$C->{$DNA_FILE}->{lbases}\t\t$C->{$DNA_FILE}->{totbases}\n\n];   
undef %{$K};

# print updated sequence file
#$DNA_FILE = qq[$outDIR/updated_sequences.fa];
#open my ($update_fh), '>>', $DNA_FILE or die "$!\n";
#for my $dna ( sort keys %DNA_SEQ )
#{  
     # all sequences marked 'u' are updated sequences that have not yet been written file 
 #   print $update_fh qq[>$dna\n$DNA_SEQ{$dna}->{seq}\n] if $DNA_SEQ{$dna}->{status} eq 'u';  #
#}
#close $update_fh;
#open update file and record sequence data summary 
#inFasta ( $DNA_FILE, \%{$K} );
#datFasta ( $DNA_FILE, \%{$K}, \%{$C} );
#print qq[$DNA_FILE\n$C->{$DNA_FILE}->{seqcount}\t\t$C->{$DNA_FILE}->{defbases}\t\t$C->{$DNA_FILE}->{nbases}\t\t$C->{$DNA_FILE}->{lbases}\t\t$C->{$DNA_FILE}->{totbases}\n\n];
undef %{$K}; 
#open scaffold file
$DNA_FILE = qq[$outDIR/new_scaffolds.fa];
inFasta ( $DNA_FILE, \%{$K} );
datFasta ( $DNA_FILE, \%{$K}, \%{$C} );
print qq[$DNA_FILE\n$C->{$DNA_FILE}->{seqcount}\t\t$C->{$DNA_FILE}->{defbases}\t\t$C->{$DNA_FILE}->{nbases}\t\t$C->{$DNA_FILE}->{lbases}\t\t$C->{$DNA_FILE}->{totbases}\n\n];   
undef %{$K};
#close $REPORT_FH;
#open my ($report), '>', qq[$outDIR/report.txt] or die "$!\n";
#print $report "$REPORT";
unlink glob "$tmpDIR/*.fa*";
rmdir "$tmpDIR";


sub scaffold_a 
{
    my $scaffold_id = shift;
    my $path_ref = shift; 
    my @scaffold;

    for my $dna ( @{$path_ref->{path}} )
    {
        if ( $dna->{strand} eq '+' )
        {  
            push @scaffold, $DNA_SEQ{$dna->{id}}->{seq};
        } 
       
        else
        {
            push @scaffold, $DNA_SEQ{$dna->{id}}->{revcom};
        }                               
              
        $DNA_SEQ{$dna->{id}}->{status} = $DNA_SEQ{$dna->{id}}->{status} eq 'o' ? 'so' : 'su';
    }
    
    my $id = $scaffold_id->();           
    my $rdna = join ( LINK, @scaffold );
    open my ($link_fh), '>>', qq[$outDIR/new_scaffolds.fa] or die "$!\n";              
    print $link_fh qq[>AX$id\n$rdna\n];
    close $link_fh;
     
}


sub scaffold_b
{
    my $scaffold_id = shift;
    my $paths_ref = shift;
    my $plus = shift;
    my $minus = shift; 
    my @scaffold;

    for my $path ( @$paths_ref )
    {
        for my $dna ( @$path )
        {
            if ( exists $plus->{$dna} )
            {   
                push @scaffold, $DNA_SEQ{$dna}->{seq};
            } 
       
            else
            {
                push @scaffold, $DNA_SEQ{$dna}->{revcom};
            }   
                                 
            $DNA_SEQ{$dna}->{status} = $DNA_SEQ{$dna}->{status} eq 'o' ? 'so' : 'su';
        }

        my $id = $scaffold_id->();           
        my $rdna = join ( LINK, @scaffold );
        open my ($link_fh), '>>', qq[$outDIR/new_scaffolds.fa] or die "$!\n";              
        print $link_fh qq[>BX$id\n$rdna\n];
        close $link_fh;
    } 
}
    


sub plus_minus 
{
    my $paths = shift; 
    my $plus = shift;
    my $minus = shift;

    for my $path ( @$paths )
    {
        if ( $path->{order} eq 'foreward' )
        {
            for my $dna ( @{$path->{path}} )
            {
                if ( $dna->{strand} eq '+' )
                {
                    $plus->{$dna->{id}}++;
                }
                
                elsif ( $dna->{strand} eq '-' )
                {
                    $minus->{$dna->{id}}++; 
                }
            }
        }
            
        elsif ( $path->{order} eq 'reverse' )
        {
            for my $dna ( @{$path->{path}} )
            {
                if ( $dna->{strand} eq '+' )
                {
                    $minus->{$dna->{id}}++;
                }
                
                elsif ( $dna->{strand} eq '-' )
                {
                    $plus->{$dna->{id}}++; 
                }
            }
        }        
    }
}


# Paths that form components should have a plausible order based on the strand positions of the vertices that they share in common   
sub order_paths 
{
    my $component_paths = shift;
    my $ordgraph = Graph::Undirected->new();
    my %vertices; 
    my %order;

    for my $path ( @$component_paths )
    {
        for ( my $i = 0; $i < scalar @{$path->{path}} - 1; $i++ )
        {
            my $v1 = @{$path->{path}}[$i]->{id};
            my $v2 = @{$path->{path}}[$i+1]->{id};
            my $s1 = @{$path->{path}}[$i]->{strand};
            my $s2 = @{$path->{path}}[$i+1]->{strand}; 
            my $f1 = "$v1".':'."$s1";
            my $r1 = $s1 eq '+' ? "$v1".':'.'-' : "$v1".':'.'+'; 
            my $f2 = "$v2".':'."$s2";
            my $r2 = $s2 eq '+' ? "$v2".':'.'-' : "$v2".':'.'+';  
            $vertices{$v1}++;
            $vertices{$v2}++;
            $ordgraph->add_edge ( $f1, $f2  );
            $ordgraph->add_edge ( $r1, $r2  );
	}
    }	 
             
    my @cc = $ordgraph->connected_components();
    return unless scalar @cc == 2;    
    my $ca = $cc[0];
    my $cb = $cc[1];  
    return unless scalar @$ca == scalar @$cb; 

    foreach ( @$ca )
    {
        my ( $dna, $strand ) = split ':', $_; 
        $order{$dna}{$strand}++;
    }        
         
    my @A = keys %order;
    my @B = keys %vertices;

    if ( scalar @A == scalar @B )    
    {
        for my $path ( @$component_paths )
        {
            foreach ( @{$path->{path}} )
            {
                if ( exists $order{$_->{id}}{$_->{strand}} )
                {
                    next;
                }

                else
                {
                    $path->{order} = 'reverse';
                    last;
                }
            }
        }
    }       
}


# remove paths whose set of vertices are included in the set of vertices in a second path 
sub set_inclusion 
{
    my $paths = shift; 
    ## Sort by length of subarray; shortest last to facilitate deletion in a loop
    @$paths = sort{ $#{ $b } <=> $#{ $a } } @$paths;
    ## Build an index to the elements in all the subarrays
    my( $i, %index ) = 0;

    for my $ref ( @$paths )
    {
        $index{ $_ } //= ++$i for @$ref;
    }
  
    ## Use the index to build bitmaps representing the subarrays
    my @bitmaps = ( '' ) x @$paths;
    for my $n ( 0 .. scalar @$paths - 1 ) 
    {
        vec( $bitmaps[ $n ], $index{ $_ }, 1 ) = 1 for @{ @$paths[ $n ] };
    }

    ## bitwise OR the bitmaps to discover wholly contained subarrays and remove them
    OUTER:for my $i ( reverse 0 .. scalar @$paths - 1 ) 
    {
        for my $j ( 0 .. $i-1 ) 
        {     
            if( ( $bitmaps[ $j ] | $bitmaps[ $i ] ) eq $bitmaps[ $j ] ) 
            {
                delete @$paths[ $i ];
                next OUTER;
            }
        }
    }

    ## remove undefs
    @$paths = grep defined, @$paths;
    return $paths; 
}


# modified depth first search to return all paths as seperate arrays from a given source vertex
# return paths are unique in that they do not include other paths as sub paths; 
# But the set of vertices in one path may include the set of vertices in a second path and so we use set operations to filter these paths ( see sub set_inclusion ) 
  
sub DFS
{
    my $start_node = shift;
    my $g = shift; 
    my @queue = ($start_node);
    my @paths;
 
    while( scalar ( @queue ) > 0 )
    {
        my $node = pop ( @queue );
        my @next_nodes;
       
        if ( index $node, ':' )
        {
            my @n = split ':', $node;
            @next_nodes = $g->successors($n[-1]);
          
            if ( ! scalar @next_nodes )
	    {
                push @paths, \@n;
                next;
            }     
        }

        else
        {
            @next_nodes = $g->successors($node);
        }
               
        @next_nodes = map { "$node".':'."$_" } @next_nodes;
        push @queue, @next_nodes;
    }
       
    return \@paths;
}



sub inFasta
{
    my $seqs_file = shift;
    my $seqs_ref = shift;
    local $/ = ">";
    open my ($sfh), '<', "$seqs_file" or die "$!\n";
    my $junk = <$sfh>;
    my $id;
    my $comment;
    while ( my $fasta = readline $sfh )
    {
        my @lines = split "\n", $fasta;
        my $def = shift @lines;
        chomp $def;
        my $seq = join ("", @lines);
        chomp $seq;
        $seq =~ s/\s//g;
        $seq =~ tr/[0-9 ]//d;
        ($id, $comment) = split (' ', $def, 2);
        $id =~ s/^\s+|\s+$//g;
        chop $id if (rindex ($id, '|') == (length $id) - 1);
        $id = substr ($id, rindex ($id, '|') + 1) if index ($id, '|');
        #$comments_ref->{$id} = $comment;
        $seqs_ref->{$id} = $seq;
    }
    close $sfh;
}

sub inPSL
{
    my $align_ref = shift;
    my $psl = shift;
    my $min_node_score = shift;
    my $min_pid = shift;
    my $dna_count = shift;
    open my ($psl_fh), '<', $psl or die "Cannot open $psl:$!\n";
    my $reg = qr /^\d+/;
    # use mean_pid to assess what scoring system to use
    # if mean pid is high and has a small range use the score with greater error weight to emphasize mismatches and qinsertions 
    my $mean_pid;
    my $mean_coverage_score; 
    my %psl_result;
    my %dd; 
    my $dna_seen = (); 

    while ( readline $psl_fh )
    {
        chomp;
        next if !/$reg/;
        my ( $matches, $mismatches, $repMatches, $ncount, $qNumInsert, $qBaseInsert, $tNumInsert, $tBaseInsert, $strand, $qName, $qSize, $qStart, $qEnd, $tName, $tSize, $tStart, $tEnd, $blockCount, $coverage_score, $high_score, $low_score, $matchTotal, $PID,  $refBlocks, $qStarts_ref, $tStarts_ref ) = parse_psl ( $_ );
        next if $PID < $min_pid;
        next if $low_score < $min_node_score;
        next if $matches < $match;               
        my $score = ( $matches + $repMatches) - ( $qNumInsert + $mismatches ); 
        
        if ( exists $psl_result{$qName}{$tName} )
	{
            my $prev_score = $psl_result{$qName}{$tName};

            if ( $prev_score < $score )
            {
                $psl_result{$qName}{$tName} = $score; 
                $align_ref->{$qName}->{$tName} = ALIGN->new ( dna => $tName, matches => $matches, mismatches => $mismatches, repmatches => $repMatches, ncount => $ncount, qinsert => $qNumInsert, qbase => $qBaseInsert, tinsert => $tNumInsert, tbase => $tBaseInsert, qsize => $qSize, tsize => $tSize, qstart => $qStart, tstart => $tStart, qend => $qEnd, tend => $tEnd, strand => $strand, bnum => $blockCount, b_sizes => $refBlocks, q_starts => $qStarts_ref, t_starts => $tStarts_ref, coverage_score => $coverage_score, high_score => $high_score, low_score => $low_score, pid => $PID );
            } 
        }

        else
	{
            $psl_result{$qName}{$tName} = $score; 
            $align_ref->{$qName}->{$tName} = ALIGN->new ( dna => $tName, matches => $matches, mismatches => $mismatches, repmatches => $repMatches, ncount => $ncount, qinsert => $qNumInsert, qbase => $qBaseInsert, tinsert => $tNumInsert, tbase => $tBaseInsert, qsize => $qSize, tsize => $tSize, qstart => $qStart, tstart => $tStart, qend => $qEnd, tend => $tEnd, strand => $strand, bnum => $blockCount, b_sizes => $refBlocks, q_starts => $qStarts_ref, t_starts => $tStarts_ref, coverage_score => $coverage_score, high_score => $high_score, low_score => $low_score, pid => $PID );   
	}
    }

    for my $rna ( keys %psl_result )
    {
        my @dna;
        push @dna, { dna => $_, score => $psl_result{$rna}{$_} } foreach keys %{$psl_result{$rna}};
        
        if ( scalar @dna > $dna_count )
        {
            @dna = sort { $b->{score} <=> $a->{score} } @dna;
            my @f = @dna[0 .. $dna_count - 1];
            $dna_seen->{$_->{dna}}++ foreach @f; 
            my %filt = map { $_->{dna} => 1 } @f; 

            for my $dna ( keys %{$align_ref->{$rna}} )
            { 
		delete $align_ref->{$rna}->{$dna} if ! exists $filt{$dna};
            }
        }  
    }


=some	    
        for ( my $i = 0; $i < $blockCount; $i++ )
        {        
            my $RNAstart = @{$align_ref->{$qName}->{$tName}->{q_starts}}[$i] + 1;
            my $DNAstart = @{$align_ref->{$qName}->{$tName}->{t_starts}}[$i] + 1;
            push @{$align_ref->{$qName}->{$tName}->{blocks}}, { id => $tName, qstart => $RNAstart, qend => $RNAstart + @{$align_ref->{$qName}->{$tName}->{b_sizes}}[$i] - 1, tstart => $DNAstart, tend => $DNAstart + @{$align_ref->{$qName}->{$tName}->{b_sizes}}[$i] - 1 };
            # why just use coverage score?  
            $align_ref->{$qName}->{$tName}->{per_cov} = percent_coverage ( \@{$align_ref->{$qName}->{$tName}->{blocks}}, $qSize, $ncount ); 
        }
    }	    
=cut
	    
    for my $rna ( keys %{$align_ref} )
    {
        my @dna = keys %{$align_ref->{$rna}};
        $dd{$_}{$rna} = () foreach @dna;
    }
                 
    my $filt = 0;
	
    for my $dna ( keys %dd )
    {
        my @rna = keys %{$dd{$dna}};
        next if scalar @rna == 1; 
        my @rel; 
        push @rel, { rna => $_, match => $align_ref->{$_}->{$dna}->{matches}, qsize => $align_ref->{$_}->{$dna}->{qsize}, tstart => $align_ref->{$_}->{$dna}->{tstart}, tend => $align_ref->{$_}->{$dna}->{tend} } foreach @rna;
        @rel = sort { $b->{match} <=> $a->{match} } @rel;
        my $top_match = $rel[0]->{match};
        my $avg_match;
        $avg_match += $_->{match} foreach @rel;
        $avg_match = $avg_match / scalar @rel;
        my @filt = grep { ( ( $_->{match} / $top_match )  < 0.10 ) && $_->{match} < 401 } @rel;
     
        for my $rna ( @filt )
        {
            delete $align_ref->{$rna->{rna}}->{$dna};
            delete $dd{$dna}{$rna->{rna}};  
            $dna_seen->{$dna}--;
            $filt++;
        }
    }
         
    for my $dna ( keys %dd )
    {
         my @rna = keys %{$dd{$dna}};
         my @scores;
         push @scores, { rna => $_, score => $psl_result{$_}{$dna} } foreach @rna;
         my @f;
         my %filt;

         if ( scalar @rna > $DEFAULT_RNA_MAX )
         {
             @scores = sort { $b->{score} <=> $a->{score} } @scores;
             @f = @scores[0..$DEFAULT_RNA_MAX - 1];
             my %filt = map { $_->{rna} => 1 } @f;
                      
            for my $rna ( @rna )
            { 
		delete $align_ref->{$rna}->{$dna} if ! exists $filt{$dna};
                $dna_seen->{$dna}--;
            }
	 }
    }
            	
    print qq[$filt dna alignments removed\n];
    my $c =()= keys %{$dna_seen}; 
    return $c;
}
 
=some        
        # insertions on the dna sequence are expected ( intronic sequence ) but insertions on query indicate poor quality alignment
        #print qq[\n$qName : $tName : $blockCount\n];
 
        #if ( $blockCount > 1 and $qNumInsert )
        #{
          # # print qq[$qNumInsert query insertions : $qBaseInsert query bases inserted\n\n]; 
           # my $qstart_prev;
         #   my $qend_prev;
            #my $tstart_prev;
            #my $tend_prev;
          
            #foreach ( @{$align_ref->{$qName}->{$tName}->{blocks}} )
           # {
                #if ( ! defined $qstart_prev )
                #{
                 #   $qstart_prev = $_->{qstart};
                  #  $qend_prev = $_->{qend};
                   # $tstart_prev = $_->{tstart};
                    #$tend_prev = $_->{tend};
                   # print qq[Query start: $_->{qstart}\tQuery end: $_->{qend}\n\n]; 
                    #next;
               # }

                #print qq[Query start: $_->{qstart}\tQuery end: $_->{qend}\n\n]; 
               # my $q_insert_start = $qend_prev + 1;
               # my $q_insert_end = $_->{qstart} - 1;
               # my $qbases = $_->{qstart} - ( $qend_prev + 1 ); 
                #print qq[$qbases qBases inserted\t$q_insert_start : $q_insert_end\n];

                # inserted bases on rna may be due to presence of divergent bases on DNA sequence.
                # it's logical to ask why there's an insertion ( ie., bases missing/added/not aligned on dna sequence, defined or gap bases? )                
                # if an insertion in alignment with rna is due to gap on dna sequence then the number of gap bases do not exactly correspond to the number of bases on the rna sequence. keep alignment only if the gap on the dna sequence varies in the number of bases from the insert on the rna by less than 10%

                #if ( $tBaseInsert and $qbases )
               # {
                    #print qq[T start: $_->{tstart}\tT end prev $tend_prev\n];
                 #   my $t_insert_start = $tend_prev + 1;
                  #  my $t_insert_end = $_->{tstart} - 1;
                   # my $tbases = $_->{tstart} - ( $tend_prev + 1 );
                    
                   # if ( $tbases )
                    #{
                        #print qq[$tbases tBases inserted\t$t_insert_start : $t_insert_end\n];
                     #   my $tseq = $strand eq '+' ? $DNA_SEQ{$tName} : rc ( $DNA_SEQ{$tName} );
                        #print qq[$tseq\n];
                        #my $l = length ($tseq);
                        #print qq[length $l\n];
                        #print qq[start $t_insert_start\tbases $tbases\n\n];
                      #  my $inserted_sequence_t = substr $tseq, $t_insert_start - 1, $tbases;  
                       # my $inserted_sequence_q = substr $RNA_SEQ{$qName}, $q_insert_start - 1, $qbases;
                        #print qq[inserted sequence q:\n$inserted_sequence_q\n]; 
                        #print qq[inserted sequence t:\n$inserted_sequence_t\n]; 
                        #my $gap_count = ( $inserted_sequence_t =~ tr/N// );

                        #if ( $gap_count == $tbases )
                        #{
                         #   my $percent_div = $qbases < $tbases ?  int $qbases / $tbases * 100 : int $tbases / $qbases * 100;
                            # insertion was due to gap on dna but gap is too large/ small for number of inserted bases
                           # if ( $percent_div >= 10 )
                          #  {
                               # print qq[query insertion due to gap : $percent_div\n]; 
                            #    delete $align_ref->{$qName}->{$tName};
                             #   last;  
                            #}

                            #else
                            #{ 
                                #print qq[query insertion due to gap : keeping $percent_div\n]; 

                            #}
                            
                        #}
                        
                        #elsif ( $tbases > 2 )
                        #{
                         #   delete $align_ref->{$qName}->{$tName};
                          #  last;
                        #}  
                   # }   

                    #print qq[$gap_count N sequence in inserted sequence t\n\n];  
                #}
 
                #elsif ( $qbases )
                #{
                    # meaning missing bases in dna sequence have broken the alignment and left unaligned rna sequence. could be due to deletion event or indicate that alignment is not biologically significant  
                    #print qq[qBase insert on query only\n\n];

                    # tolerate small deletions in coding sequence ( insertions in rna sequence ) 
                 #   if ( $qbases > 2 )
                   # {
                  #      delete $align_ref->{$qName}->{$tName};
                    #    last;
                    #}
                #}

               # $qstart_prev = $_->{qstart};
                #$qend_prev = $_->{qend};
		#$tstart_prev = $_->{tstart};
		#$tend_prev = $_->{tend};                        
           # }
        #}
        # any gap bases included in alignment represent none covered rna sequence; inserted bases are by definition outside block start/end positions 

=cut

sub rc
{
    my $seq = shift; 
   ( my $rc = reverse $seq ) =~ tr /ATCGatcg/TAGCtagc/;
    return $rc;
}


sub percent_coverage
{
    my $ss_ref = shift;
    my $rna_length = shift;
    my $ncount = shift;

    @$ss_ref = sort { $a->{qstart} <=> $b->{qstart} || $a->{qend} <=> $b->{qend} } @$ss_ref;
    my $start = 1;
    my $x = 0;

    foreach ( @$ss_ref )
    {
        $x += $_->{qstart} - $start if $start < $_->{qstart};
        $start = $_->{qend} if $_->{qend} > $start;
    }

    $x += $rna_length - $start if $start != $rna_length;
    my $percent_missing_sequence = ( $x - $ncount ) / $rna_length * 100;
    my $c = 100 - $percent_missing_sequence;
    $c = sprintf ("%.1f", $c );
    return $c;
}


sub parse_psl
{
    my ( @b, @t, @q );
    my ($line) = @_;
    my @result = split(/\t/,$line);
    my $matches = $result[0];     # Number of bases that matches the query matched to the target
    my $misMatches = $result[1];  # Number of bases that don't match
    my $repMatches = $result[2];  # Number of bases that match but are part of repeats
    my $nCount = $result[3];      # Number of 'N' bases
    my $qNumInsert = $result[4];  # Number of inserts in query
    my $qBaseInsert = $result[5]; # Number of bases inserted in query
    my $tNumInsert = $result[6];  # Number of inserts in target
    my $tBaseInsert = $result[7]; # Number of bases inserted in target
    my $strand = $result[8];      # '+' or '-' for query strand. For translated alignments, second '+'or '-' is for genomic strand
    my $qName = $result[9];       # Query sequence name
    my $qSize = $result[10];      # Query sequence size
    my $qStart = $result[11];     # Alignment start position in query
    my $qEnd = $result[12];       # Alignment end position in query
    my $tName = $result[13];      # Target sequence name
    my $tSize = $result[14];      # Target sequence size
    my $tStart = $result[15];     # Alignment start position in target
    my $tEnd = $result[16];       # Alignment end position in target
    my $blockCount = $result[17]; # Number of blocks in the alignment (a block contains no gaps)
    my $blockSizes = $result[18]; # Comma-separated list of sizes of each block
    my $qStarts = $result[19];    # Comma-separated list of starting positions of each block in query
    my $tStarts = $result[20];    # Comma-separated list of starting positions of each block in target
    
    return if ! $blockSizes; 
   
    $blockSizes =~ s/[, ]$//;
    $tStarts =~ s/[, ]$//;
    @b = split /,/, $blockSizes; #
    @t = split /,/, $tStarts;
    @q = split /,/, $qStarts;
    #$qName =~ s/^\s+|\s+$//g;
    #$sName =~ s/^\s+|\s+$//g;
 
    foreach ( $qName, $tName )
    {
        chop $_ if (rindex ($_, '|') == (length $_) - 1);
        $_ = substr ($_, rindex ($_, '|') + 1) if index ($_, '|');
    }

    
    my $p = ( ( $misMatches + $qNumInsert ) / $matches ) * 100;  
    my $score = ( ( $matches + $misMatches + $repMatches ) / $qSize ) * 100;
    $score = sprintf ( "%.1f", $score );  

    my $match_total  = $matches + $misMatches + $repMatches;
    my $percent_id = sprintf ("%.1f", ( ( ( $matches + $repMatches ) / $match_total ) * 100 ) );
    my $coverage_score = ( $match_total / $qSize ) * 100;
  
   
    my $mm = $misMatches * 2; 
    my $penalty = ( ( $misMatches + $qNumInsert ) / ( $matches - $mm - $qNumInsert ) ) * 100;
    my $alignment_score_high =  sprintf ( "%.1f", $coverage_score - $penalty );

    $mm = $misMatches * 3;
    $penalty =  ( ( $misMatches + $qNumInsert ) / ( $matches - $mm - $qNumInsert ) ) * 100;
    my $alignment_score_low = sprintf ( "%.1f", $coverage_score - $penalty );  
    $coverage_score = sprintf ("%.1f", $coverage_score);
    
    # even if strand is minus we keep block coordinates as plus

    if ( $strand eq '-' )
    {
        @b = reverse @b;
        my @reverse_t = reverse @t;
        my @reverse_q = reverse @q;       

        for ( my $i = 0; $i < $blockCount; $i++ )
        {
            my $tstart = $tSize - ( $b[$i] + $reverse_t[$i] );     
            $t[$i] = $tstart; 
            my $qstart = $qSize - ( $b[$i] + $reverse_q[$i] );
            $q[$i] = $qstart;
        } 
    }  
    
    return( $matches,$misMatches,$repMatches,$nCount,$qNumInsert,$qBaseInsert,$tNumInsert,$tBaseInsert,$strand,$qName,$qSize,$qStart,$qEnd,$tName,$tSize,$tStart,$tEnd,$blockCount, $score, $alignment_score_high, $alignment_score_low, $match_total, $percent_id, \@b, \@q, \@t );
}

sub groupSequences 
{
    my $align_ref = shift;
    my $g = Graph::Undirected->new();
    
    for my $rna ( keys %{$align_ref} ) 
    {        
        for my $dna ( keys %{$align_ref->{$rna}} ) 
        {
	    $g->add_edge("rna:$rna", "dna:$dna");
        }
    }

    my @groups;

    for my $raw_group ( $g->connected_components() ) 
    {
        my %group = ( dna => [], rna => [] );

        for ( @$raw_group ) 
        {
	    my ($type, $val) = split(/:/, $_, 2);
	    push @{ $group{$type} }, $val;
        }

        push @groups, \%group; 
    }
   
    return \@groups;
}


sub RDassemble
{
    my $dna_ref = shift; 
    my $seen_ref = shift;
    my $seq_id = shift; 
    my $DNA;
    $DNA->{$_} = () foreach @$dna_ref;
     
    for my $q ( @$dna_ref )
    {       
        next unless exists $DNA->{$q};        
      
	for my $s ( @$dna_ref )
	{
            next if exists $seen_ref->{$q}->{$s} || exists $seen_ref->{$s}->{$q} || $q eq $s;
            $seen_ref->{$q}->{$s} = (); 
            next if ( ( $DNA_SEQ{$q}->{gap} >= 10 and $DNA_SEQ{$s}->{gap} >= 10 ) || ( ! $DNA_SEQ{$q}->{gap} and ! $DNA_SEQ{$s}->{gap} >=10 ) );    
            next unless exists $DNA->{$s};
            my @hsp;
            my ( $strand, $c ) = BL2SEQ ( $q, $s, \@hsp );
            next if ! $c || $c > 2;   
            my $s_seq;
            
            if ( $strand eq '+-' )
            {
                foreach ( @hsp )
                {
                    $_->{sstart} = $DNA_SEQ{$s}->{length} - $_->{sstart} + 1;  
                    $_->{send} = $DNA_SEQ{$s}->{length} - $_->{send} + 1;  
                }

                $s_seq = $DNA_SEQ{$s}->{revcom}; 
            } 
            
            $s_seq = $DNA_SEQ{$s}->{seq} unless $s_seq;                         
            my $new_seq = ();
            my $k;

            if ( $c == 1 )  
	    { 
                my $align_q = $hsp[0]->{qend} - ( $hsp[0]->{qstart} - 1 );  
                my $align_s = $hsp[0]->{send} - ( $hsp[0]->{sstart} - 1 );
                my $len_q = $DNA_SEQ{$q}->{length};
                my $len_s = $DNA_SEQ{$s}->{length};
                my $red; 

                if ( $hsp[0]->{percentID} >= 99 ) 
                {
                    $red = ( $len_q == $align_q ) ? $q : ( $len_s == $align_s ) ? $s : '';  
                         
                    if ( $red ) 
                    {
			delete $ALIGN{$_}{$red} foreach @{$DNA_SEQ{$red}->{rna}};                    
                        $DNA_SEQ{$red}->{status} = ( $DNA_SEQ{$red}->{status} eq 'o' ) ? 'or' : 'ur';
                        delete $DNA->{$red}; 
                        print $REPORT_FH qq[Red\t$q\t$s\t\t\t$red\n];
                        last if $red eq $q;   
			next if $red eq $s;
                    }
                }
 
                elsif ( $hsp[0]->{qstart} - 1 == ( $len_q - $align_q ) && $hsp[0]->{sstart} == 1 ) 
                {
                    next; 
                    #unless $align_q > 90;
                    #my $c = substr ( $s_seq, $align_s ); 
                    #$c = $DNA_SEQ{$q}->{seq}.$c;
                    #$new_seq = \$c;     
                }

                elsif ( $hsp[0]->{sstart} - 1 == ( $len_s - $align_s ) && $hsp[0]->{qstart} == 1 )          
                {
                    next;
                    # unless $align_s > 90;
                    #my $c = substr ( $DNA_SEQ{$q}->{seq}, $align_q ); 
                    #$c = $s_seq.$c;
                    #$new_seq = \$c; 
                }
            }
            
            ( $k, $new_seq ) = merge2seq ( \@hsp ) if ! defined $$new_seq; 
             
            if ( defined $new_seq )
            {
                my $id = $seq_id->();
                $DNA->{$id} = ();  
                ( my $rc = reverse $$new_seq ) =~ tr /ATCGatcg/TAGCtagc/; 
                $DNA_SEQ{$id} = DNA->new ( seq => $$new_seq, revcom => $rc, status => 'u', length => length $$new_seq, gap => $$new_seq =~ tr/N//, file => 1 );
                newSequenceFile ( $new_seq, $tmpDIR, "$id", '>' );
                print $REPORT_FH qq[End Overlap\t$q\t$s\t$id\n] if ! $k;
                print $REPORT_FH qq[Gap Close\t$q\t$s\t$id\t$k\n] if $k;
               
                $k = ();
                # for each removed dna obtain non redundant list of RNA and delete from queue 
                for my $dna ( $q, $s )
                { 
                    delete $DNA->{$dna};  
                    push @{$DNA_SEQ{$id}->{predecessors}}, $_ foreach @{$DNA_SEQ{$dna}->{predecessors}};
                   
                    for my $rna ( @{$DNA_SEQ{$dna}->{rna}} )
                    { 
                        $k->{$rna} = (); 
                        delete $ALIGN{$rna}{$dna};
                    }

                    $DNA_SEQ{$dna}->{status} = ( $DNA_SEQ{$dna}->{status} eq 'o' ) ? 'ou' : 'uu';
                }
                # set expected new seq rna alignments         
                foreach ( keys %{$k} )
                {
                    $ALIGN{$_}{$id} = 'update'; 
                    push @{$DNA_SEQ{$id}->{rna}}, $_;
	        }
                
                my @dna = keys %{$DNA}; 
		return RDassemble ( \@dna, $seen_ref, $seq_id );
            }
        }   
    }
}



sub BL2SEQ
{
    my ( $q, $s, $hsp_ref ) = @_;
    my $hsp_count = 0;
    my $blast = `/work/lianglab/bin/blast-2.2.26/bin/bl2seq -i"$tmpDIR/$q.fa" -j"$tmpDIR/$s.fa" -pblastn -D1 -e1e-20 -FF`;
    open my ($blast_fh), '<', \$blast;
    my $filter_strand; 
    

  LINE:
    while (<$blast_fh>)
    {
        next if (index ($_, '#') == 0);
        chomp;
        my ( $percent_identity, $align_length, $mismatches, $qstart, $qend, $sstart, $send, $strand, $gap_openings, $bit_score, $evalue ) = parseBlast ( $_ );

        if ( ! $hsp_count && ( $percent_identity < 99 || $align_length < 30 ) )
        {
            return 0;
        }

        elsif ( ! $hsp_count && ( $percent_identity >= 99 && $align_length >= 30 ) )
        {
            push @$hsp_ref, { query => $q, subject => $s, str => $strand,  qstart => $qstart, qend => $qend, sstart => $sstart, send => $send, alignlength => $align_length, mismatches => $mismatches, percentID => $percent_identity, gap => $gap_openings, bit => $bit_score, evalue => $evalue };
            $hsp_count++;
            $filter_strand = $strand; 
            next;
        }

        elsif ( $hsp_count && ( $percent_identity >= 98 && $strand eq $filter_strand ) )
        {
            for my $block ( @$hsp_ref )
            {
                next LINE unless ( ( $qstart > $block->{qend} || $qend < $block->{qstart} ) &&  ( $sstart > $block->{send} || $send < $block->{sstart} ) ); 
            }

            push @$hsp_ref, { query => $q, subject => $s, str => $strand, qstart => $qstart, qend => $qend, sstart => $sstart, send => $send, alignlength => $align_length, mismatches => $mismatches, percentID => $percent_identity, gap => $gap_openings, bit => $bit_score, evalue => $evalue };
            $hsp_count++;
        }
    }
    close $blast_fh;
    return ( $filter_strand, $hsp_count );
}


sub parseBlast
{
    my $line              = shift;
    my @result            = split(/\s+/,$line);
    my $query             = $result[0];
    my $subject           = $result[1];
    my $percent_identity  = $result[2];
    my $align_length      = $result[3];
    my $mismatches        = $result[4];
    my $gap_openings      = $result[5];
    my $query_start       = $result[6]; 
    my $query_end         = $result[7];
    my $subject_start     = $result[8];
    my $subject_end       = $result[9];
    my $e_value           = $result[10];
    my $bit_score         = $result[11]; 
    my $strand            = ( $subject_start < $subject_end ) ? '++' : '+-'; 
    return ($percent_identity, $align_length, $mismatches, $query_start, $query_end, $subject_start, $subject_end, $strand, $gap_openings, $bit_score, $e_value);  
}


sub merge2seq
{
    my $hsp_ref = shift; 
    my $count = scalar @$hsp_ref; 
    my $strand = @$hsp_ref[0]->{str};  
    @$hsp_ref   = sort { $a->{qstart} <=> $b->{qstart} } @$hsp_ref;
    my $qstart  = @$hsp_ref[0]->{qstart};
    my $qid     = @$hsp_ref[0]->{query};     
    my $qend    = @$hsp_ref[$count-1]->{qend};
    @$hsp_ref   = sort { $a->{sstart} <=> $b->{sstart} } @$hsp_ref;
    my $sstart  = @$hsp_ref[0]->{sstart};
    my $sid     = @$hsp_ref[0]->{subject}; 
    my $send    = @$hsp_ref[$count-1]->{send};
    my $qseq = $DNA_SEQ{$qid}->{seq};
    my $sseq;
    $sseq = $DNA_SEQ{$sid}->{seq} if $strand eq '++';
    $sseq = $DNA_SEQ{$sid}->{revcom} if $strand eq '+-';
    my @new;   
    my $mismatch = 0;
    my $Nfill = 0;
    my $s;
    my $e;
    my $end; 
    my $start = 0; 
 
    if ( $qstart > $sstart )
    {
        $start = $qstart - $sstart; 
        $s = substr $qseq, 0, $start, ''; 
        
        for my $hsp ( @$hsp_ref )
        { 
            $hsp->{qstart} = $hsp->{qstart} - $start;
            $hsp->{qend} = $hsp->{qend} - $start;
        } 
    }

    elsif ( $sstart > $qstart )
    {
	$start = $sstart - $qstart; 
        $s = substr $sseq, 0, $start, '';
   
        for my $hsp ( @$hsp_ref )
        { 
            $hsp->{sstart} = $hsp->{sstart} - $start;
            $hsp->{send} = $hsp->{send} - $start;
        }
    }
   
    if ( $DNA_SEQ{$sid}->{length} - $send > $DNA_SEQ{$qid}->{length} - $qend )
    {
        $end = ( ( $DNA_SEQ{$sid}->{length} - $send ) - ( $DNA_SEQ{$qid}->{length} - $qend ) );
        $e = substr $sseq, -$end, $end, '';  
    }
    
    elsif ( $DNA_SEQ{$qid}->{length} - $qend > $DNA_SEQ{$sid}->{length} - $send )
    {
        $end = ( ( $DNA_SEQ{$qid}->{length} - $qend ) - ( $DNA_SEQ{$sid}->{length} - $send ) );
        $e = substr $qseq, -$end, $end, '';
    }   
 
    my @subject = split '', $sseq;
    my @query   = split '', $qseq;
    
   
    my $i = 0;
    my $j = 0;
    my $c = 0;
    
    for my $hsp ( @$hsp_ref )
    { 
        my $qstart = $hsp->{qstart};
        my $qend   = $hsp->{qend};
        my $sstart = $hsp->{sstart};
        my $send   = $hsp->{send}; 
        my $gap = $hsp->{gap}; 
        my $n = 0;
        while ( ( $i < $qstart - 1 ) and ( $j < $sstart - 1 ) ) 
        {
            if ( $query[$i] eq $subject[$j] )
            {
                if ( $subject[$j] eq 'N' )
                {
                    $new[$c] = 'N';
                }

                else 
                {
                    $new[$c] = $subject[$j];
                }
            }

            elsif ( $subject[$j] ne $query[$i] )
            {
                if  ( $subject[$j] eq 'N' && $query[$i] ne 'N'  )
                {
                    $new[$c] = $query[$i]; 
                    $Nfill++;
                }
          
                elsif ( $subject[$j] ne 'N' && $query[$i] eq 'N' ) 
		{
		    $new[$c] = $subject[$j]; 
                    $Nfill++;
                }
                    
                else
                { 
                    $new[$c] = $subject[$j]; 
		    $mismatch++;
                }
            }
            return 0 if $mismatch; 
            $i++; $j++; $c++;
	}
        
        while ( $i < $qstart - 1 )
        {
            $n++ if $query[$i] ne 'N'; 
            $i++;
            return 0 if $n;  
        }
         
        while ( $j < $sstart - 1 ) 
        {
            $n++ if $subject[$j] ne 'N';
            $j++;
            return 0 if $n;  
            
        }
 
        while ( $i < $qend )   
        {
            $new[$c] = $query[$i]; 
            $i++; $c++; 
	}

        while ( $j < $send )
        {
            $j++;
        }
    }
    
    while ( defined $subject[$j] && defined $query[$i] )
    {
        if ( $subject[$j] eq $query[$i] )
        {
            if ( $subject[$j] eq 'N' )
            {
                $new[$c] = 'N';
            }

            else 
            {
                $new[$c] = $subject[$j];
            }
        }

        elsif ( $subject[$j] ne $query[$i] )
        {
            if  ( $subject[$j] eq 'N' and $query[$i] ne 'N'  )
            { 
                $new[$c] = $query[$i]; 
                $Nfill++;
            }
          
            elsif ( $subject[$j] ne 'N' and $query[$i] eq 'N' ) 
	    {
		$new[$c] = $subject[$j]; 
                $Nfill++;
            }
                    
            else
            {
                $new[$c] = $subject[$j]; 
		$mismatch++;
            }
	}
        return 0 if $mismatch;  
        $i++; $j++; $c++; 
    }
      
    if ( ! $mismatch and $Nfill )
    {
        my $new_sequence = join '', @new;
        $new_sequence = $s.$new_sequence if $s;
        $new_sequence = $new_sequence.$e if $e;
	return ( $Nfill, \$new_sequence );   
    }

    return 0;
}


sub newSequenceFile
{
    my $sequence_ref = shift; 
    my $folder_name = shift; 
    my $file_name = shift;
    my $mode = shift; 
    open my ($new_seq), $mode, "$folder_name/$file_name.fa" or die "$!\n";
    print $new_seq ">$file_name\n$$sequence_ref"; 
    close $new_seq; 
}


sub idGenerate 
{
    my $header = shift;
    my $group_id = shift;  
    my $c = 0;

    return sub {
	$c++;
	my $name = $header.'_'.$group_id.'_'.$c;   
	return $name;
    } 
}


sub formatFasta 
{
    my $seq_ref = shift;
    my $len = 80;
    my $formatted_seq;

    while ( my $chunk = substr( $$seq_ref, 0, $len, "" ) )
    {
        $formatted_seq .= "$chunk\n"; 
    }
    
    return \$formatted_seq; 
}



# takes a fasta file name and corresponding sequence data and generates assembly metrics  
sub datFasta
{
    my $file_name = shift;
    my $seq_ref = shift;
    my $dat_ref = shift;
    #$file_name = substr ( $file, 0, ( index $file, '.fa') ) if ( index $file, '.fa' ); 
    my $count =()= keys %{$seq_ref};
    $dat_ref->{$file_name}->{seqcount} = $count;
    $count = 0;
    # N bases
    my $ncount = 0;
    # linker bases
    my $lcount = 0;
    my $dcount = 0;
    while ( my ( $id, $sequence ) = each %{$seq_ref} )
    {
        $count += length $sequence;
        $dcount += $sequence =~ tr/ATCG//; 
        $ncount += $sequence =~ tr/N//;
        $lcount += $sequence =~ tr/n//;
    }
     
    $dat_ref->{$file_name}->{totbases} = $count; 
    $dat_ref->{$file_name}->{defbases} = $dcount;
    $dat_ref->{$file_name}->{nbases} = $ncount;
    $dat_ref->{$file_name}->{lbases} = $lcount;  
}
