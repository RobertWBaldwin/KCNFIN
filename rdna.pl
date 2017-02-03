#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use Graph; 
use RDseq;
use Carp; 
use Getopt::Long;
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


my ( $DNA_FILE, $RNA_FILE, $PSL_FILE, $PATH_SCORE_MIN, $NODE_SCORE_MIN, $PID_MIN, $DNA_MAX, $GAP_PID_MAX, $GAP_DNA_MAX, $INTRON_MAX, $SCONNECT, $OVERLAP, $scaffold, $gap, $redundancy, $help ); 
# default values for scaffold options
my $DEFAULT_TMP_DIR;
my $DEFAULT_OUT_DIR; 
my $DEFAULT_PATH_SCORE = 50;
my $DEFAULT_NODE_SCORE = 4.1;
my $DEFAULT_PID        = 98.1;
my $DEFAULT_DNA_MAX    = 200;
my $DEFAULT_RNA_MAX    = 201; 
my $DEFAULT_INTRON_MAX = 100000;
# default values for gap options
my $DEFAULT_GAP_PID    = 97.6; 
my $DEFAULT_GAP_DNA    = 15;
my $DEFAULT_OVERLAP    = 8;
my $DEFAULT_SCONNECT   = 400;

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
    "sconnect=i" => \$SCONNECT,
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
$SCONNECT = defined $SCONNECT ? $SCONNECT : $DEFAULT_SCONNECT;
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

my $zz =()= keys %DNA_SEQ;
print qq[$zz DNA SEQUENCES MARKED O\n];


my %Z;
my %Y;
inFasta ( $DNA_FILE, \%Z );
my %O;

for my $s ( keys %Z )
{
    if ( exists $DNA_SEQ{$s} and $DNA_SEQ{$s}->{status} eq 'o' )
    {
        $O{$s} = $Z{$s};
        delete $Z{$s};
    }
}

my $ccount  = 0;
my $ddcount = 0;
my $nncount = 0;
my $llcount = 0;

while ( my ( $id, $sequence ) = each %Z )
{
        $ccount += length $sequence;
        $ddcount += $sequence =~ tr/ATCG//; 
        $nncount += $sequence =~ tr/N//;
        $llcount += $sequence =~ tr/n//;
}

my $f =()= keys %Z;
print qq[$f original assembly sequences exlcuded\n];
print qq[Defined_Bases: $ddcount\nN_bases:$nncount\n\n\n];

$f =()= keys %O;
$ccount  = 0;
$ddcount = 0;
$nncount = 0;
$llcount = 0;

while ( my ( $id, $sequence ) = each %O )
{
        $ccount += length $sequence;
        $ddcount += $sequence =~ tr/ATCG//; 
        $nncount += $sequence =~ tr/N//;
        $llcount += $sequence =~ tr/n//;
}

print qq[$f original sequences included\n];
print qq[Defined_Bases: $ddcount\nN_bases: $nncount\n\n\n];
# get stats for all 'o' sequences as base line 
# get stats for all absent sequences as base line 

$K = ();
$C = ();
@groups = sort { @{$a->{rna}} <=> @{$b->{rna}} || @{$a->{dna}} <=> @{$b->{dna}} } @groups; 

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
        my $v = scalar @dna; 
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
                next if grep { $_->{qstart} >= $alignment_a->{qend} - $OVERLAP && $_->{qend} <= $alignment_b->{qstart} + $OVERLAP } @{$RNA_SEQ{$rna}->{alignments}}; 
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
        my $E = $g->edges();        
       # print qq[Edges: $E\n]; 
        next if $E > 125;
        # determine all possible scaffolding paths       
        for my $source_node ( @source )
        {
           my $dfs_paths = DFS( $source_node, $g );
           push @paths, $_ foreach @$dfs_paths; 
        }

        
        my $p = scalar @paths;
        #print qq[Paths: $p\n];             
        # remove any paths whose set of vertices is present in a longer path 
        my $filt_paths = ( scalar @paths > 1 ) ? set_inclusion(\@paths) : \@paths;  
        my @path_score;
        #print qq[Scored paths\n\n];        
        
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
        #my $max_path_b = exists $path_score[1] ? $path_score[1] : 0;
        @{$RNA_SEQ{$rna}->{alignments}} = sort { $b->{low_score} <=> $a->{low_score} } @{$RNA_SEQ{$rna}->{alignments}};
        # skip RNA whose max scoring path is too low
        if ( $max_path_a->{score} < $PATH_SCORE_MIN )
        {
            $RNA_SEQ{$rna}->{status} = 's';
            next;
        }
        # skip RNA whose max scoring path is less than the top scoring alignment + constant     
        next if $max_path_a->{score} < ${$RNA_SEQ{$rna}->{alignments}}[0]->{low_score} + 3;
        #next if defined $max_path_b->{score} and $max_path_a->{score} == $max_path_b->{score}; 
        my $path = PATH->new ( rna => $rna, score => $max_path_a->{score},  );
        my %strands = map { $_->{dna} => { str => $_->{strand}, scr => $_->{low_score} } } @{$RNA_SEQ{$rna}->{alignments}};
        push @{$path->{path}}, { id => $_, strand => $strands{$_}->{str}, score => $strands{$_}->{scr} } foreach @{$max_path_a->{path}};
        push @MAX_PATHS, $path;  
    }

    @E = ();
    @RNA = ();
    my $c = scalar @MAX_PATHS;
    next unless $c;
    my %plus;
    my %minus;
  
    # group has one max path   
    if ( scalar @MAX_PATHS == 1 )
    {
        plus_minus(\@MAX_PATHS, \%plus, \%minus);

        for my $path ( @MAX_PATHS )
        {
            my @scaffold;
            my $id = $scaffold_id->(); 
            push @scaffold, $_->{id} foreach @{$path->{path}};
            one ( \@scaffold, \%plus, $i );
            scaffold_paths(\@scaffold, \%plus, $id) if scalar @scaffold;    
        }

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
    
            plus_minus(\@component_paths, \%plus, \%minus);

            for my $path ( @component_paths )
            {
                my @scaffold;
                my $id = $scaffold_id->(); 
                push @scaffold, $_->{id} foreach @{$path->{path}};
                one ( \@scaffold, \%plus, $i );
                scaffold_paths(\@scaffold, \%plus, $id) if scalar @scaffold;
            }

            next;
        }
    
                        
        # components here are "complex" consisting of multiple max paths 
        # first step is to find the most plausible order for these paths based on the strand positions of the dna sequences that they have in common (+/-)
        plus_minus ( \@component_paths, \%plus, \%minus );
        my @pm = grep { exists $plus{$_} && exists $minus{$_} } keys %vertices;          

        if ( scalar @pm )
        {
            # If you pass this condition then some of the paths will have been reversed  
            order_paths ( \@component_paths );
            %plus = ();
            %minus = ();            
            plus_minus ( \@component_paths, \%plus, \%minus ); 
            @pm = grep { exists $plus{$_} && exists $minus{$_} } keys %vertices;
            next if scalar @pm; # you need a plausible set of path orders to proceed   
        }
       
        # and now actually reverse the paths as per the new strand orientations       
        for my $path ( @component_paths )
        {
            if ( $path->{order} eq 'reverse' )
            {
                my @reverse = reverse @{$path->{path}};
                @{$path->{path}} = ();
                push @{$path->{path}}, $_ foreach @reverse;
                $path->{order} = 'foreward'; 
            }
        }            

        my $digraph = Graph->new(); 

        for my $path ( @component_paths )
        {
            for ( my $i = 0; $i < scalar @{$path->{path}} - 1; $i++ )
            {
                if ( ! $digraph->has_edge( @{$path->{path}}[$i]->{id}, @{$path->{path}}[$i+1]->{id} ) )
		{
                    $digraph->add_edge ( @{$path->{path}}[$i]->{id}, @{$path->{path}}[$i+1]->{id} );
		}
            }
        } 
   
        next if $digraph->is_cyclic;
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
            for my $scaffold ( @$filt_paths )
            {
                my $id = $scaffold_id->();
                one ( $scaffold, \%plus, $i );
                scaffold_paths($scaffold, \%plus, $id) if scalar @$scaffold;
            }

            next;
        }

        ##### No consensus path was found #######

        $K = ();
        @raw_paths = ();
        @source = (); 
        $C = ();        
            
        for my $path ( @component_paths )
        {
            
            my %match = map { $_->{dna} => $_->{matches} } @{$RNA_SEQ{$path->{rna}}->{alignments}}; 
            for ( my $i = 0; $i < scalar @{$path->{path}} - 1; $i++ )
            {
                my $a = @{$path->{path}}[$i]->{id};
		my $b = @{$path->{path}}[$i+1]->{id};
                $K->{$a}->{$b} = (); 
                push @raw_paths, { edge => [$a, $b], match => [ $match{$a}, $match{$b} ] };  
	    }
        }
 
        $K = ();
        @raw_paths = sort { @{$b->{match}}[0] <=> @{$a->{match}}[0] || @{$b->{match}}[1] <=> @{$a->{match}}[1]  } @raw_paths;

        foreach ( @raw_paths )
        {
            my $a = @{$_->{edge}}[0];
            my $b = @{$_->{edge}}[1];
            my $ma = @{$_->{match}}[0];
            my $mb = @{$_->{match}}[1];
            next if exists $C->{$a}; 
            $C->{$a}->{$b} = ();
        } 

        @raw_paths = sort { @{$b->{match}}[1] <=> @{$a->{match}}[1] || @{$b->{match}}[0] <=> @{$a->{match}}[0]  } @raw_paths;
        
        foreach ( @raw_paths )
        {
            my $a = @{$_->{edge}}[0];
            my $b = @{$_->{edge}}[1];
            my $ma = @{$_->{match}}[0];
            my $mb = @{$_->{match}}[1];
            next if exists $K->{$b};  
            $K->{$b}->{$a} = [$mb,$ma];   
        }
        
        my $final_graph = Graph->new();
        @source = ();
        @raw_paths = ();

        for my $a ( keys %{$C} )
        {
            for my $b ( keys %{$C->{$a}} )
            {
                $final_graph->add_edge ( $a, $b ) if exists $K->{$b}->{$a} and @{$K->{$b}->{$a}}[0] > 100 and @{$K->{$b}->{$a}}[1] > 100; 
	    }
        } 
	         
        next if $final_graph->is_cyclic;
        # recompute the paths
        @source = $final_graph->source_vertices();

        for my $source_vertex ( @source )
        {
            my $dfs = DFS ( $source_vertex, $final_graph );
            push @raw_paths, $_ foreach @$dfs; 
        }
       
        my $final_paths = ( scalar @raw_paths > 1 ) ? set_inclusion( \@raw_paths ) : \@raw_paths;
        @$final_paths = sort { $a->[0] cmp $b->[0] || $a->[1] cmp $b->[1] || $a->[2] cmp $b->[2] || $a->[3] cmp $b->[3] } @$final_paths;
     
        for my $scaffold ( @$final_paths )
        {
            my $id = $scaffold_id->(); 
            one ( $scaffold, \%plus, $i );
            scaffold_paths($scaffold, \%plus, $id) if scalar @$scaffold;
        }
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
$K = ();

# load remainder dna file
inFasta ( $DNA_FILE, \%{$K} );
# record sequence data summary
datFasta ( $DNA_FILE, \%{$K}, \%{$C} );
print qq[$DNA_FILE\n$C->{$DNA_FILE}->{seqcount}\t\t$C->{$DNA_FILE}->{defbases}\t\t$C->{$DNA_FILE}->{nbases}\t\t$C->{$DNA_FILE}->{lbases}\t\t$C->{$DNA_FILE}->{totbases}\n\n];  
 
$K = ();

$DNA_FILE = qq[$outDIR/updated_sequences.fa];
open my ($update_fh), '>', $DNA_FILE or die "$!\n";


for my $dna ( sort keys %DNA_SEQ )
{
    # any sequence still marked 'u' was not scaffolded   
    print $update_fh qq[>$dna\n$DNA_SEQ{$dna}->{seq}\n] if $DNA_SEQ{$dna}->{status} eq 'u';
}
close $update_fh;
#open update file and record sequence data summary 
inFasta ( $DNA_FILE, \%{$K} );
datFasta ( $DNA_FILE, \%{$K}, \%{$C} );
print qq[$DNA_FILE\n$C->{$DNA_FILE}->{seqcount}\t\t$C->{$DNA_FILE}->{defbases}\t\t$C->{$DNA_FILE}->{nbases}\t\t$C->{$DNA_FILE}->{lbases}\t\t$C->{$DNA_FILE}->{totbases}\n\n];
$K = (); 
#open scaffold file
$DNA_FILE = qq[$outDIR/new_scaffolds.fa];
inFasta ( $DNA_FILE, \%{$K} );
datFasta ( $DNA_FILE, \%{$K}, \%{$C} );
print qq[$DNA_FILE\n$C->{$DNA_FILE}->{seqcount}\t\t$C->{$DNA_FILE}->{defbases}\t\t$C->{$DNA_FILE}->{nbases}\t\t$C->{$DNA_FILE}->{lbases}\t\t$C->{$DNA_FILE}->{totbases}\n\n];   


$K = ();
$C = ();
$ccount  = 0;
$ddcount = 0;
$nncount = 0;
$llcount = 0;
$f = ();
my %Org;
my %Upd;
my %UU;
my %UO;
my %SO;
my %SU;
$f =()= keys %DNA_SEQ;
print qq[$f total sequences in DNA_SEQ\n];

for my $s ( keys %DNA_SEQ )
{
    if ( $DNA_SEQ{$s}->{status} eq 'o' )
    {
        $Org{$s} = $DNA_SEQ{$s}->{seq};
        delete $DNA_SEQ{$s};
    }    

    elsif ( $DNA_SEQ{$s}->{status} eq 'u' )
    {
        $Upd{$s} = $DNA_SEQ{$s}->{seq};
        delete $DNA_SEQ{$s};
    }

    elsif ( $DNA_SEQ{$s}->{status} eq 'uo' )
    {
        $UO{$s} = $DNA_SEQ{$s}->{seq};
        delete $DNA_SEQ{$s};
    }
    
    elsif ( $DNA_SEQ{$s}->{status} eq 'uu' )
    {
        $UU{$s} = $DNA_SEQ{$s}->{seq};
        delete $DNA_SEQ{$s};
    }

    elsif ( $DNA_SEQ{$s}->{status} eq 'so' )
    {
        $SO{$s} = $DNA_SEQ{$s}->{seq};
        delete $DNA_SEQ{$s};
    }

    elsif ( $DNA_SEQ{$s}->{status} eq 'su' )
    {
        $SU{$s} = $DNA_SEQ{$s}->{seq};
        delete $DNA_SEQ{$s};
    }

    else
    {
        print qq[$s has not status\n];
    }

}

while ( my ( $id, $sequence ) = each %Org )
{
        $ccount += length $sequence;
        $ddcount += $sequence =~ tr/ATCG//; 
        $nncount += $sequence =~ tr/N//;
        $llcount += $sequence =~ tr/n//;
}

$f =()= keys %Org;
print qq[$f DNA sequences Org\nDef_Bases: $ddcount\nN_bases: $nncount\n\n\n];

$ccount  = 0;
$ddcount = 0;
$nncount = 0;
$llcount = 0;
$f =()= keys %Upd;
while ( my ( $id, $sequence ) = each %Upd )
{
        $ccount += length $sequence;
        $ddcount += $sequence =~ tr/ATCG//; 
        $nncount += $sequence =~ tr/N//;
        $llcount += $sequence =~ tr/n//;
}

print qq[$f DNA sequences Upd\nDef_Bases: $ddcount\nN_bases: $nncount\n\n\n];



$ccount  = 0;
$ddcount = 0;
$nncount = 0;
$llcount = 0;
$f =()= keys %SO;
while ( my ( $id, $sequence ) = each %SO )
{
        $ccount += length $sequence;
        $ddcount += $sequence =~ tr/ATCG//; 
        $nncount += $sequence =~ tr/N//;
        $llcount += $sequence =~ tr/n//;
}

print qq[$f DNA sequences SO\nDef_Bases: $ddcount\nN_bases: $nncount\n\n\n];





$ccount  = 0;
$ddcount = 0;
$nncount = 0;
$llcount = 0;
$f =()= keys %SU;
while ( my ( $id, $sequence ) = each %SU )
{
        $ccount += length $sequence;
        $ddcount += $sequence =~ tr/ATCG//; 
        $nncount += $sequence =~ tr/N//;
        $llcount += $sequence =~ tr/n//;
}

print qq[$f DNA sequences SU\nDef_Bases: $ddcount\nN_bases: $nncount\n\n\n];




$ccount  = 0;
$ddcount = 0;
$nncount = 0;
$llcount = 0;
$f =()= keys %UO;
while ( my ( $id, $sequence ) = each %UO )
{
        $ccount += length $sequence;
        $ddcount += $sequence =~ tr/ATCG//; 
        $nncount += $sequence =~ tr/N//;
        $llcount += $sequence =~ tr/n//;
}

print qq[$f DNA sequences UO\nDef_Bases: $ddcount\nN_bases: $nncount\n\n\n];





$ccount  = 0;
$ddcount = 0;
$nncount = 0;
$llcount = 0;
$f =()= keys %UU;
while ( my ( $id, $sequence ) = each %UU )
{
        $ccount += length $sequence;
        $ddcount += $sequence =~ tr/ATCG//; 
        $nncount += $sequence =~ tr/N//;
        $llcount += $sequence =~ tr/n//;
}

print qq[$f DNA sequences UU\nDef_Bases: $ddcount\nN_bases: $nncount\n\n\n];


close $REPORT_FH;
open my ($report), '>', qq[$outDIR/report.txt] or die "$!\n";
#print $report "$REPORT";
unlink glob "$tmpDIR/*.fa*";
rmdir "$tmpDIR";



sub scaffold_paths
{
    my $path = shift;
    my $plus = shift;
    my $id = shift;
    my @scaffold;
 
    for my $dna ( @$path )
    {
        if ( exists $plus->{$dna} || $DNA_SEQ{$dna}->{status} eq 'u' )
        {   
            push @scaffold, $DNA_SEQ{$dna}->{seq};
        } 

        else
        {
            push @scaffold, $DNA_SEQ{$dna}->{revcom};
        }   

        $DNA_SEQ{$dna}->{status} = $DNA_SEQ{$dna}->{status} eq 'o' ? 'so' : 'su';
    }
    
    my $rdna = join ( LINK, @scaffold );
    open my ($link_fh), '>>', qq[$outDIR/new_scaffolds.fa] or die "$!\n";              
    print $link_fh qq[>$id\n$rdna\n];
    close $link_fh;     
}


sub one
{
    my $scaffold = shift;
    my $plus = shift;
    my $i = shift;
    my $id = idGenerate ( 'RDNAU', $i );
    my $connect = Graph->new();
    my $k;
        
    for ( my $i = 0; $i < scalar @$scaffold - 1; $i++ )
    {
        $connect->add_edge(@$scaffold[$i],@$scaffold[$i+1]);
    }

    for my $dna ( @$scaffold ) 
    {
        open my ( $os_fh ), '>', qq[$tmpDIR/$dna.fa] or die "$!\n";

        if ( $plus->{$dna} )
        {
            print $os_fh qq[>$dna\n$DNA_SEQ{$dna}->{seq}\n];
            $DNA_SEQ{$dna}->{file} = '+';
        }
    
        else
        {
            print $os_fh qq[>$dna\n$DNA_SEQ{$dna}->{revcom}\n];
            $DNA_SEQ{$dna}->{file} = '-';
        }   

        close $os_fh; 
    }
    
    check_connections( $connect, $id, $k );
    my @edges = $connect->edges();
    # no updated sequences; scaffold is the same
    if ( scalar @edges == scalar @$scaffold - 1 )
    {
        return;
    }
    # no edges then no scaffolding to do; will remain 'u' and be printed to updated sequence file
    elsif ( ! scalar @edges )
    {
        @$scaffold = ();
        return;
    }
   
    elsif ( scalar @edges == 1 )
    {
        @$scaffold = ();
        my $edge = shift @edges;
        @$scaffold[0] = @$edge[0];
        @$scaffold[1] = @$edge[1];
    }
    
    else
    {
        @$scaffold = (); 
        my @source = $connect->source_vertices();
        my $path = DFS ( $source[0], $connect );
        my $path_ref = shift @$path;
        push @$scaffold, $_ foreach @$path_ref;
    }

}    


sub check_connections
{
    my $connect = shift; 
    my $seq_id = shift; 
    my $count = shift;
    my @edges = $connect->edges();
    my $new_seq;
    my $k;
    my $c;
    
    for my $edge ( @edges )
    {
        my $v1 = @$edge[0];
        my $v2 = @$edge[1];
        next if exists $count->{$v1}->{$v2}; 
        $count->{$v1}->{$v2} = ();
 
        my @hsp;
        my ( $strand, $c ) = BL2SEQ ( $v1, $v2, \@hsp );
        next unless $c;           
        next if $strand eq '+-'; 
        $new_seq = $c == 1 ? overlap(\@hsp) : undef;
        ( $k, $new_seq ) =  merge2seq(\@hsp) if ! defined $$new_seq && ( $DNA_SEQ{$v1}->{gap} || $DNA_SEQ{$v2}->{gap} );
        
        if ( defined $$new_seq )
        {
            my $id = $seq_id->();
            $id = ! exists $DNA_SEQ{$id} ? $id : '';

            if ( ! $id )
            { 
                while ( $id = $seq_id->() )
                {
                    last if ! exists $DNA_SEQ{$id};
                }
            }  

            ( my $rc = reverse $$new_seq ) =~ tr /ATCGatcg/TAGCtagc/; 
            $DNA_SEQ{$id} = DNA->new ( seq => $$new_seq, status => 'u', length => length $$new_seq, gap => $$new_seq =~ tr/N//, file => 1 );
            my $status = $DNA_SEQ{$id}->{status};
            newSequenceFile ( $new_seq, $tmpDIR, "$id", '>' );
            my $align_q = $hsp[0]->{qend} - ( $hsp[0]->{qstart} - 1 );
            my @remove = grep { @$_[0] eq $v1 || @$_[1] eq $v1 || @$_[0] eq $v2 || @$_[1] eq $v2 } @edges;            

            for my $dna ( $v1, $v2 )
            { 
                $DNA_SEQ{$dna}->{status} = $DNA_SEQ{$dna}->{status} eq 'o' ? 'uo' : 'uu';
            }
            
            for my $edge ( @remove ) 
            {
                $connect->delete_edge(@$edge[0], @$edge[1]);

                if ( @$edge[0] eq $v2 and @$edge[1] ne $v1 )
                {
                    $connect->add_edge($id, @$edge[1]);
                }

                elsif ( @$edge[1] eq $v1 and @$edge[0] ne $v2 )
                {
                    $connect->add_edge(@$edge[0], $id);
                }
            }
            
            return if ! $connect->edges();
            return check_connections ($connect, $seq_id, $count );
        }
    }
}
        

sub overlap 
{
    my $hsp = shift; 
    my $query = @$hsp[0]->{query};
    my $subject = @$hsp[0]->{subject};
    my $align_q = @$hsp[0]->{qend} - ( @$hsp[0]->{qstart} - 1 );  
    my $align_s = @$hsp[0]->{send} - ( @$hsp[0]->{sstart} - 1 );
    my $len_q = $DNA_SEQ{$query}->{length};
    my $len_s = $DNA_SEQ{$subject}->{length};
    my $query_sequence   = $DNA_SEQ{$query}->{file} eq '+' ? $DNA_SEQ{$query}->{seq} : $DNA_SEQ{$query}->{revcom};
    my $subject_sequence = $DNA_SEQ{$subject}->{file} eq '+' ? $DNA_SEQ{$subject}->{seq} : $DNA_SEQ{$subject}->{revcom};    
    my $new_sequence; 
       
    if ( @$hsp[0]->{qstart} - 1 == ( $len_q - $align_q ) && @$hsp[0]->{sstart} == 1 ) 
    {
        my $x = substr ( $query_sequence, $len_q - $align_q - 200 );
        my $y = substr ( $subject_sequence, $align_s, 200 );
        my $xy = $x.$y;
        print qq[END OVERLAP:\n$query <---> $subject $align_q base pairs\nregion:\n$xy\n\n];
        my $c = substr ( $subject_sequence, $align_s );
        $c = $query_sequence.$c;
        $new_sequence = \$c;     
    }

    return $new_sequence
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
    #my $rna_count = shift; 
    open my ($psl_fh), '<', $psl or die "Cannot open $psl:$!\n";
    my $reg = qr /^\d+/;
    # use mean_pid to assess what scoring system to use
    # if mean pid is high and has a small range use the score with greater error weight to emphasize mismatches and qinsertions 
    my $mean_pid;
    my $mean_coverage_score; 
    my %psl_result;
    my %dd; 
    my $dna_seen = (); 
    my %hundred; 

    while ( readline $psl_fh )
    {
        chomp;
        next if !/$reg/;
        my ( $matches, $mismatches, $repMatches, $ncount, $qNumInsert, $qBaseInsert, $tNumInsert, $tBaseInsert, $strand, $qName, $qSize, $qStart, $qEnd, $tName, $tSize, $tStart, $tEnd, $blockCount, $coverage_score, $high_score, $low_score, $matchTotal, $PID,  $refBlocks, $qStarts_ref, $tStarts_ref ) = parse_psl ( $_ );
        next if $PID < $min_pid;
        next if exists $hundred{$qName}; 
        next if $low_score < $min_node_score;
        next if $tSize < 100 || $matchTotal < 50; 
        $hundred{$qName} = () if $low_score >= 95;   
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
        my @filt = grep { ( ( $_->{match} / $top_match )  < 0.15 ) && $_->{match} < 400 } @rel;
     
        for my $rna ( @filt )
        {
            delete $align_ref->{$rna->{rna}}->{$dna};
            delete $dd{$dna}{$rna->{rna}};  
            $dna_seen->{$dna}--;
            $filt++;
        }
    }
     
    %dd = ();
 
    for my $rna ( keys %{$align_ref} )
    {
        my @dna = keys %{$align_ref->{$rna}};
        $dd{$_}{$rna} = () foreach @dna;
    }
    
    for my $dna ( keys %dd )
    {
         my @rna = keys %{$dd{$dna}};
         my $length = length($DNA_SEQ{$dna}); 
         my @scores;
         push @scores, { rna => $_, score => $psl_result{$_}{$dna} } foreach @rna;
         my @f;
         my %filt;
         my $rna_max = $length <= 1000 ? 50 : $length <= 10000 ? 100 : $length <= 50000 ? 200 : 300;   

         if ( scalar @rna > $rna_max )
         {
             @scores = sort { $b->{score} <=> $a->{score} } @scores;
             @f = @scores[0..$rna_max - 1];
             my %filt = map { $_->{rna} => 1 } @f;
                      
            for my $rna ( @rna )
            { 
		delete $align_ref->{$rna}->{$dna} if ! exists $filt{$dna};
                $dna_seen->{$dna}--;
            }
	 }
    }
            	
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

=some
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
=cut


sub BL2SEQ
{
    my ( $q, $s, $hsp_ref ) = @_;
    my $hsp_count = 0;
    my $blast = `/work/lianglab/bin/blast-2.2.26/bin/bl2seq -i"$tmpDIR/$q.fa" -j"$tmpDIR/$s.fa" -pblastn -D1 -e1e-20`;
    open my ($blast_fh), '<', \$blast;
    my $filter_strand; 
    
  LINE:
    while ( readline $blast_fh )
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
    my $count   = scalar @$hsp_ref; 
    my $query   = @$hsp_ref[0]->{query};
    my $query_sequence = $DNA_SEQ{$query}->{file} eq '+' ? $DNA_SEQ{$query}->{seq} : $DNA_SEQ{$query}->{revcom};  
    my $subject = @$hsp_ref[0]->{subject};
    my $subject_sequence = $DNA_SEQ{$subject}->{file} eq '+' ? $DNA_SEQ{$subject}->{seq} : $DNA_SEQ{$subject}->{revcom};
    
    @$hsp_ref   = sort { $a->{qstart} <=> $b->{qstart} } @$hsp_ref;
    my $qstart  = @$hsp_ref[0]->{qstart};
    my $qend    = @$hsp_ref[$count-1]->{qend};
      
    @$hsp_ref   = sort { $a->{sstart} <=> $b->{sstart} } @$hsp_ref;
    my $sstart  = @$hsp_ref[0]->{sstart};
    my $send    = @$hsp_ref[$count-1]->{send};
    
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
        $s = substr $query_sequence, 0, $start, ''; 
        
        for my $hsp ( @$hsp_ref )
        { 
            $hsp->{qstart} = $hsp->{qstart} - $start;
            $hsp->{qend} = $hsp->{qend} - $start;
        } 
    }

    elsif ( $sstart > $qstart )
    {
	$start = $sstart - $qstart; 
        $s = substr $subject_sequence, 0, $start, '';
   
        for my $hsp ( @$hsp_ref )
        { 
            $hsp->{sstart} = $hsp->{sstart} - $start;
            $hsp->{send} = $hsp->{send} - $start;
        }
    }
   
    if ( $DNA_SEQ{$subject}->{length} - $send > $DNA_SEQ{$query}->{length} - $qend )
    {
        $end = ( ( $DNA_SEQ{$subject}->{length} - $send ) - ( $DNA_SEQ{$query}->{length} - $qend ) );
        $e = substr $subject_sequence, -$end, $end, '';  
    }
    
    elsif ( $DNA_SEQ{$query}->{length} - $qend > $DNA_SEQ{$subject}->{length} - $send )
    {
        $end = ( ( $DNA_SEQ{$query}->{length} - $qend ) - ( $DNA_SEQ{$subject}->{length} - $send ) );
        $e = substr $query_sequence, -$end, $end, '';
    }   
 
    my @subject = split '', $subject_sequence;
    my @query   = split '', $query_sequence;
    
   
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
        my $q = $DNA_SEQ{$query}->{file} eq '+' ? $DNA_SEQ{$query}->{seq} : $DNA_SEQ{$query}->{revcom};
        my $s = $DNA_SEQ{$subject}->{file} eq '+' ? $DNA_SEQ{$subject}->{seq} : $DNA_SEQ{$subject}->{revcom};
        print qq[GAP CLOSED $Nfill N bases:\n>$query\n$q\n\n>$subject\n$s\n\n>new_sequence\n$new_sequence\n\n];  
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
