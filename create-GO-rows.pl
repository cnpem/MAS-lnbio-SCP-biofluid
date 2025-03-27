#!/usr/bin/perl -w
# Copyright 2025 Guilherme P. Telles.
# This program is distributed under the terms of WTFPL v2.

use Getopt::Long;
use LWP::Simple;
use LWP::UserAgent;
use Time::HiRes('usleep');
use XML::LibXML;
use XML::LibXML::XPathContext;

$| = 1;
$" = ', ';

my $tabf = '';
my $ecof = '';
my $qcolumn = 1;
my $gof = '';
my $datad = '.';
my $function = 'count';
my $level = -1;
my $allev = 0;
my $help = 0;

my $verb = 1;

Getopt::Long::Configure('bundling');

my $st = GetOptions('i|input=s' => \$tabf,
		    'b|eco=s' => \$ecof,
		    'q|column=i' => \$qcolumn, 
		    'l|level=i' => \$level,		 
		    'g|go=s' => \$gof,
		    'f|function=s' => \$function,
		    'a|all_evidence' => \$allev,
		    'w|data=s' => \$datad,
		    'h|help' => \$help);

#print STDERR "i => $tabf, b => $ecof, qcolumn => $qcolumn, l => $level, g => $gof,",
#"f => $function, a => $allev, w => $datad\n";
#print STDERR "st => $st, @ARGV\n";

if (!$st || @ARGV != 0 || $tabf eq '' || $ecof eq '' || $help) {
  
  $help = 1;
}


if ($help) {
  print<<_END_;

Create rows of attributes based on GO annotations.

The input is a tsv file with protein quantifications in rows and
conditions in columns.

The input file must have protein ids in the first column (separated by
semicolon).  Optionally, the second and other consecutive columns may
have metadata.  The columns after the last metadata column must contain
quantifications.
  
The output to stdio is a tsv with the aggregation of Biological
Process GOs in rows and conditions in columns (the same as input),
having the GO code in the first column.

If the aggregation function is count, then each entry
[GO,condition] in the output is the number of proteins that were
quantified in the condition and annotated with that GO in Uniprot.

If the quantification function is sum, then each entry [GO,condition]
in the output is the sum of quantifications of proteins that are
present in the condition and annotated with that GO in Uniprot.

usage: create-GO-rows.pl -i file -b file 

arguments:
  -i  A tsv file with proteins in rows and conditions in columns.
  -b  An Evidence Ontology file in OBO format.

options:
  -q  The number of columns preceding the first column with
      quantification data.  Defaults to 1.
  -l  If provided, the ancestor level to which each GO will be
      generalized, 0 meaning the root of GO hierarchy.
  -g  A tree GO Ontology file in OBO format.  Mandatory if -l is given.
  -f  The aggregation function for GO.  May be count or sum.
      Default is count.
  -a  Use all GO annotations in Uniprot.  Default is to use only GO
      annotations with experimental evidence.
  -w  A directory where Uniprot files (one for each protein) will be
      written to and read from.  The program will try to download
      Uniprot files that are not in that directory.
      If not given, the current directory will be used.

_END_

  exit(2);
}

(-d $datad) or die("$datad: $!.");

($ecoref,$econamesref) = obo_load_dag($ecof);

if ($level >= 0) {
  ($goref,$gonamesref) = obo_load_tree($gof);
}

my $TAB;
open($TAB,"<",$tabf) or die("$tabf: $!");

# Header:
$_ = <$TAB>;
chomp;
my @header = split(/\t/);

my $ncol = @header;

my %T = ();  # The term for each GO.
my %G = ();  # An array of counts for each GO.

my %ancestor = ();

while (<$TAB>) {

  chomp;
  my @F = split(/\t/);
  
  my @acc = split(/;/,"$F[0]",2);
  
  my $prot = $acc[0];

  if ($prot eq '') {
    next;
  }
  
  $verb and print STDERR "$prot: ";

    
  ### Download UniProt xml or use local files:
  $down = uprot_download_xml($prot,$datad);
  
  if ($down eq 'exists') {
    $verb and print STDERR "already stored locally\n";
  }
  elsif ($down eq 'fail') {
    $verb and print STDERR "unable to download from UniProt\n";
    next;
  }
  elsif ($down ne $prot) {
    $verb and print STDERR "downloaded from UniProt as $down\n";
  }
  else {
    $verb and print STDERR "downloaded from UniProt\n";
  }
  

  ### Get GO:
  $unif = "$datad/$prot.xml";
  @go = ();
  @go = uprot_get_go($unif,$ecoref,$econamesref);

  for $h (@go) {

    %H = %{$h};
    
    $go = $H{GO};
    $term = $H{term};

    #print STDERR "$prot, $H{GO}, $H{term}, $H{evidence_name}, $H{evidence_is_experimental}\n";

    if ($H{term} =~ /^P:/ && ($allev || $H{evidence_is_experimental} eq 'yes')) {  
      
      if ($level >= 0) {

	if (exists($ancestor{$go})) {
	  # Use ancestor GO searched previously:
	  $go = $ancestor{$go};
	  $term = %{$gonamesref}{$go};
	}
	else {
	  # Search in GO-tree:
	  my @path = ();
	  $st = obo_path_tree($goref,$go,'GO:0008150',\@path);

	  if ($st) {
	    if (@path > $level) {
	      $ancestor{$go} = $path[@path-1-$level];
	      $go = $path[@path-1-$level];
	      $term = %{$gonamesref}{$go};
	    }
	    else {
	      $ancestor{$go} = $go;
	    }
	  }
	  else {
	    # Failed to find the ancestor, disregard this GO.
	    next;
	  }
	}
      }      

      if (!exists(${T}{$go})) {
	${T}{$go} = $term;
	$a = $ncol-1;
	$G{$go} = [ (0) x $a ];
      }
	
      for ($i=$qcolumn; $i<$ncol; $i++) {
	
	if ($F[$i] eq 'NA') {
	  $F[$i] = 0;
	}
	
	if ($F[$i] > 0) {
	  if ($function eq 'count') {
	    $G{$go}[$i-1] += 1;
	  }
	  elsif ($function eq 'sum') {
	    $G{$go}[$i-1] += $F[$i];
	  }
	  else {
	    die("Unknown GO aggregation function\n");
	  }
	}
      }
    }
  }
}

close($TAB);

$" = "	";

for $go (keys %G) {
  print "$go $T{$go}\t@{$G{$go}}\n";
}


$verb and print STDERR "done.\n";



################################################################################
# (hash-ref,hash-ref) obo_load_dag($obo-file)
#
# Load an ontology from a file in OBO format.
#
# Return
# a reference to a hash that contains the dag as term-id => @is_a and
# a reference to a hash that maps term-id => $name.

sub obo_load_dag {

  my $filename = shift;

  my $OBO;
  open($OBO,"<",$filename) || die("Unable to open $filename");

  my %dag = ();  # Every node has an array of parent ids.
  my @name = ();
  
  my $id = '';
  my @isa = ();
  my $name = '';
  my $in = 0;
  my $obsolete = 0;

  while (<$OBO>) {
    chomp;

    /^\[/ && do {

      if ($in) {
	$dag{$id} = [ @isa ];
	$name{$id} = $name;
      }
      
      $id = '';
      @isa = ();
      $name = '';
      $obsolete = 0;
      
      $in = /^\[Term\]\s*$/ ? 1 : 0;
      next;
    };

    $in && /^id:\s*([0-9A-Za-z:]+)\s*/ && do {
      $id = $1;
      next;
    };


    $in && /^is_obsolete: true\s*$/ && do {
      $obsolete = 1;
      $name = "$name : obsolete";
      next;
    };

    $in && $obsolete && /^replaced_by:\s*([0-9A-Za-z:]+)\s*/ && do {
      $name = "$name : $id replaced by $1";
      next;
    };

    
    $in && /^name: (.*)$/ && do {
      $name = $1;
    };
    
    $in && /^is_a:\s*([0-9A-Za-z:]+)\s*/ && do {
      push(@isa,$1);
      next;
    }
  }

  if ($in) {
    $dag{$id} = [ @isa ];
    $name{$id} = $name;
  }
  
  close($OBO);
  
  return (\%dag,\%name);
}



################################################################################
# (hash-ref,hash-ref) obo_load_tree($obo-file)
#
# Load an OBO file with an ontology that is a tree.
#
# Obsolete terms are loaded with name augmented with the prefix
# ': obsolete : id replaced by id' if a replaced_by exists.
#
# Return
# a reference to a hash that contains the tree as term-id => $is_a and
# a reference to a hash that maps term-id => $name.

sub obo_load_tree {

  my $filename = shift;

  my $OBO;
  open($OBO,"<",$filename) || die("Unable to open $filename");

  my %tree = ();  # Every node stores its single parent id.
  my @name = ();
  
  my $id = '';
  my $name = '';
  my $isa = '';
  my $in = 0;
  my $obsolete = 0;

  while (<$OBO>) {
    chomp;

    /^\[/ && do {

      if ($in) {
	$tree{$id} = $isa;
	$name{$id} = $name;
      }
      
      $id = '';
      $isa = '';
      $name = '';
      $obsolete = 0;

      $in = /^\[Term\]\s*$/ ? 1 : 0;
      next;
    };

    $in && /^id:\s*([0-9A-Za-z:]+)\s*/ && do {
      $id = $1;
      next;
    };

    $in && /^is_obsolete: true\s*$/ && do {
      $obsolete = 1;
      $name = "$name : obsolete";
      next;
    };

    $in && $obsolete && /^replaced_by:\s*([0-9A-Za-z:]+)\s*/ && do {
      $name = "$name : $id replaced by $1";
      next;
    };
    
    $in && /^name: (.*)$/ && do {
      $name = $1;
      next;
    };
    
    $in && /^is_a:\s*([0-9A-Za-z:]+)\s*/ && do {
      $isa = $1;
      next;
    }
  }

  if ($in) {
    $tree{$id} = $isa;
    $name{$id} = $name;
  }
  
  close($OBO);
  
  return (\%tree,\%name);
}



################################################################################
# obo_descends_dag(\%dag, $id, $ancestor)
#
# Return 1 if id "descends" from ancestor, and return 0 otherwise.

sub obo_descends_dag {

  my $dagref = shift;
  my $curr = shift;
  my $from = shift;

  if ($curr eq $from) {
    return 1;
  }
  
  for my $child (@{%{$dagref}{$curr}}) {
    if ($child eq $from || obo_descends_dag($dagref,$child,$from)) {
      return 1;
    }
  }

  return 0;
}



################################################################################
# obo_path_tree(\%tree, $id, $ancestor, \@path)
#
# Enumerate the nodes in the path from id to ancestor, appending to @path.
# Return 1 if a path exists, 0 otherwise.

sub obo_path_tree {

  my $treeref = shift;
  my $curr = shift;
  my $from = shift;
  my $pathref = shift;

  my %tree = %{$treeref};

  while ($curr && $curr ne $from) {
    push(@{$pathref},$curr);
    $curr = $tree{$curr};
  }

  if ($curr) {
    push(@{$pathref},$curr);
    return 1;
  }
  else {
    while (@{$pathref}) {
      pop(@{$pathref});
    }
    return 0;
  }
}



################################################################################
# string uprot_download_xml($accession,$directory)
#
# Download a UniProt entry by accession as xml.
# If it fails, it tries to download a synonym entry.
# If it fails again, it creates a placeholder file.
#
#
# If $directory/$accession.xml exists already then it returns "exists".
#
# If the entry itself was downloaded then it saves the entry as
# $directory/$accession.xml and returns $accession.
#
# If a synonym $directory/X.xml exists already or was downloaded then
# it creates a symbolic link from $directory/$accession.xml to
# $directory/X.xml, and returns X.
#
# If the entry or a synonym fail to download then it creates an empty
# file $directory/$accession-fail.xml and returns "fail".
#
# The directory defaults to the cwd.  It dies on file operations failure.

sub uprot_download_xml {

  my $acc = shift;
  my $dir = shift;

  !defined($dir) && ($dir = '.');

  my $get = $acc;  
  my $retry = 0;
  
  while (1) {

    if (-f "$dir/$get.xml") {
      if ($retry) {
	if (!-f "$dir/$acc.xml") {
	  symlink("$dir/$get.xml","$dir/$acc.xml") or die("$!\n");
	}
	return $get;
      }      
      return "exists";	 
    }
    
    # Try downloading xml:
    my $url = "https://www.ebi.ac.uk/proteins/api/proteins/$get";
    my $agent = LWP::UserAgent->new(agent => 'libwww-perl gpt@ic.unicamp.br');
    my $response = $agent->get("$url", 'Accept'=>'application/xml');
    usleep(10000);

    while (my $wait = $response->header('Retry-After')) {
      sleep($wait);
      $response = $agent->get($response->base);
    }
  
    if ($response->is_success) {

      open(my $fh,">","$dir/$get.xml") || die("$! : $dir/$get.xml\n");
      print $fh $response->content."\n";
      close($fh);
      
      if ($retry) {
	if (!-f "$dir/$acc.xml") {
	  symlink("$dir/$get.xml","$dir/$acc.xml") or die("$!\n");
	}
	return $get;
      }

      return $acc;
    }
    else {
      
      if (!$retry) {
	# Try getting a synonym acc. On success, retry a download. On failure, give up:    
	my @syn = uprot_download_synonyms($acc);
	$retry = 1;

	my $i = 0;
	while ($i < @syn && $syn[$i] eq $acc) {
	  $i++;
	}
	
	if ($i < @syn) {
	  $get = $syn[$i];
	  next;
	}
      }

      # Failed at getting a synonym or downloading using a synonym acc:
      touch("$dir/$acc-fail.xml");
      return 'fail';
    }
  }
}



################################################################################
# array-of-hashes uprot_get_go($xml-file, \%eco-dag, \%eco-names) 
#
# Get GO annotations from Uniprot XML.
#
# Return an array of hashes with fields type, description, at,
# evidence, experimental (0 or 1).
  
sub uprot_get_go {

  my $filename = shift;
  my $dagr = shift;
  my $namesr = shift;

  my %dag = %{$dagr};
  my %names = %{$namesr};  

  my @G = ();

  # Process xml:
  my $dom = XML::LibXML->load_xml(location => $filename);
  my $xpc = XML::LibXML::XPathContext->new($dom);
  $xpc->registerNs('up','https://uniprot.org/uniprot');

  
  # Collect data on evidences in parallel arrays:
  foreach my $ev ($xpc->findnodes('//up:entry/up:dbReference')) {

    my %h = ();

    my $type = $ev->findnodes('./@type');
    my $key = $ev->findnodes('./@id');

    if ($type eq "GO") {
      $h{$type} = $key;

      foreach my $prop ($xpc->findnodes('./up:property',$ev)) {

	my $type = $prop->findnodes('./@type');
	my $value = $prop->findnodes('./@value');
	$h{$type} = $value;

	if ($type eq "evidence") {
	  $h{evidence_name} = $names{$value};
	  if (obo_descends_dag(\%dag,$value,'ECO:0000006')) {
	    $h{evidence_is_experimental} = 'yes';
	  }
	  else {
	    $h{evidence_is_experimental} = 'no';
	  }	    
	}
      }

      # print "=>$_ $h{$_}\n" for (keys %h);
      push(@G, { %h } );    
    }
  }

  return @G;
}
