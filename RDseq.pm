#!/usr/bin/perl -w

package PATH;
use Moose;
has 'order', is => 'rw', default => 'foreward', isa => 'Str';
has 'rna', is => 'rw', isa => 'Str';
has 'path', is => 'rw', isa => 'ArrayRef[Node]';
has 'score', is => 'rw', isa => 'Num';

package Node;
use Moose;
has 'id', is => 'rw', isa => 'Str';
has 'strand', is => 'rw', isa => 'Str';
has 'score', is => 'rw', isa => 'Num';


package RDseq;
use Moose;
has 'seq', is => 'rw', isa => 'Str'; 
has 'length', is => 'rw', isa => 'Int'; 
has 'status', is => 'rw', isa => 'Str';

package BLOCKS;
use Moose;
has 'id', is => 'rw', isa => 'Str'; 
has 'dna', is => 'rw', isa => 'Str';
has 'qstart', is => 'rw', isa => 'Int';
has 'qend', is => 'rw', isa => 'Int';
has 'sstart', is => 'rw', isa => 'Int';
has 'send', is => 'rw', isa => 'Int'; 
has 'length', is =>'rw', isa => 'int';

package ALIGN;
use Moose; 
has 'rna', is => 'rw', isa => 'Str';
has 'dna', is => 'rw', isa => 'Str';
has 'matches', is => 'rw', isa => 'Int';
has 'mismatches', is => 'rw', isa => 'Int';
has 'repmatches', is => 'rw', isa => 'Int';
has 'coverage_score', is => 'rw', isa => 'Num';
has 'high_score', is => 'rw', isa => 'Num';
has 'low_score', is => 'rw', isa => 'Num';
has 'pid', is => 'rw', isa => 'Num';
has 'ncount', is => 'rw', isa => 'Int';
has 'qinsert', is => 'rw', isa => 'Int';
has 'qbase', is => 'rw', isa => 'Int';
has 'qsize', is => 'rw', isa => 'Int';
has 'tinsert', is => 'rw', isa => 'Int';
has 'tbase', is => 'rw', isa => 'Int';
has 'tsize', is => 'rw', isa => 'Int'; 
has 'qstart', is => 'rw', isa =>'Int';
has 'qend', is => 'rw', isa => 'Int';
has 'tstart', is => 'rw', isa => 'Int';
has 'tend', is => 'rw', isa => 'Int';
has 'strand', is => 'rw', isa => 'Str';
has 'bnum', is => 'rw', isa => 'Int';
has 'b_sizes', is => 'rw', isa => 'ArrayRef[Int]';
has 'q_starts', is => 'rw', isa => 'ArrayRef[Int]';
has 't_starts', is => 'rw', isa => 'ArrayRef[Int]';
has 'pcov', is => 'rw', isa => 'Num';
has 'blocks', is => 'rw', isa => 'ArrayRef[BLOCKS]'; 


package DNA;     
use Moose;
extends 'RDseq';
has 'revcom', is => 'rw', isa => 'Str';
has 'gap',    is => 'rw', isa => 'Int';
has 'rna',    is => 'rw', isa => 'ArrayRef[Str]';
has 'predecessors', is => 'rw', isa => 'ArrayRef[Str]'; 
has 'file', is => 'rw', default => '0', isa => 'Bool'; 
has 'scaffolded', is => 'rw', default =>'0', isa => 'Bool';


package RNA;
use Moose;
extends 'RDseq';
has 'dna', is => 'rw', isa => 'ArrayRef[BLOCKS]';
has 'coverage', is => 'rw', isa => 'Int';
has 'alignments', is=> 'rw', isa => 'ArrayRef[ALIGN]';

1;

