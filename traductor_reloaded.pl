#!perl
use Getopt::Long;
my %opts = ();
GetOptions(\%opts,'i=s', 'n=i', 'o=s', 'h');
if(($opts{'h'})||(scalar(keys(%opts))==0)){
print <<HELP;
NAME
        conectador.pl

PARAMETTERS 

	

AUTHOR
        Roberto Galindo Ramirez
        LCG 4a Generacion
	Instituto de Fisiologia Celular
	Ciudad Universitaria UNAM

DESCRIPTION
        


EXAMPLE
     perl traductor.pl -i fragmentoorfgef_alignment.aln -n 29 
 
	
        
FORMAT 
        


HELP
}
else{
####GENETIC CODE
%tabla = (
        'TAA'=>'.','TAG'=>'.','TGA'=>'.',
        'GCT'=>'A','GCC'=>'A','GCA'=>'A','GCG'=>'A',
        'TGT'=>'C','TGC'=>'C',
        'GAT'=>'D','GAC'=>'D',
        'GAA'=>'E','GAG'=>'E',
        'TTT'=>'F','TTC'=>'F',
        'GGT'=>'G','GGC'=>'G','GGA'=>'G','GGG'=>'G',
        'CAT'=>'H','CAC'=>'H',
        'ATT'=>'I','ATC'=>'I','ATA'=>'I',
        'AAA'=>'K','AAG'=>'K',
        'CTT'=>'L','CTC'=>'L','CTA'=>'L','CTG'=>'L',
        'TTA'=>'L','TTG'=>'L',
        'ATG'=>'M',
        'AAT'=>'N','AAC'=>'N',
        'CCT'=>'P','CCC'=>'P','CCA'=>'P','CCG'=>'P',
        'CAA'=>'Q','CAG'=>'Q',
        'CGT'=>'R','CGC'=>'R','CGA'=>'R','CGG'=>'R',
        'AGA'=>'R','AGG'=>'R',
        'TCT'=>'S','TCC'=>'S','TCA'=>'S','TCG'=>'S',
        'AGT'=>'S','AGC'=>'S',
        'ACT'=>'T','ACC'=>'T','ACA'=>'T','ACG'=>'T',
        'GTT'=>'V','GTC'=>'V','GTA'=>'V','GTG'=>'V',
        'TGG'=>'W',
        'TAT'=>'Y','TAC'=>'Y',
);

####
	$file=$opts{'i'};
	$nseq=$opts{'n'};
	open(ARCHI,$file)or die"\nFailure to open\n";
	@alin=<ARCHI>;
	close(ARCHI);
#####PARTE DE CONCATENADO DE SECUENCIAS
	for($k=3;$k<=4+$nseq;$k++){
		$seqs[$k-3]=substr($alin[$k],0,21);
	}
	for($i=0;$i<scalar(@alin);$i++){
		if($alin[$i]=~/\*/){
			if($alin[$i-1]=~/ATG/){
				$index=$i;
				last;	
			}
		}
	}
	for($i=$index;$i<scalar(@alin);$i+=$nseq+2){
		for($j=$i-$nseq,$k=0;$j<$i;$j++,$k++){
			$seqs[$k].=substr($alin[$j],22);
		}
	}
#####PARTE DE CONCATENADO DE SECUENCIAS


	for($k=0,$proteina[$k]='';$k<$nseq;$k++){
		for($i=22;$i<length($seqs[$k])-2;$i++){
			$metionina=substr($seqs[$k],$i,3);
			if($metionina eq "ATG"){
				for($j=$i;$j<length($seqs[$k]);$j+=3){
					$codon=substr($seqs[$k],$j,3);
					$aa=$tabla{$codon};
					if($aa ne "."){
					$proteina[$k].=$aa;
					}
					else{
						last;
					}
				}
				#print"$proteina[$k]\n";
				$comienzo[$k]=$i;
				last;
			}
		}
	}
####La parte donde se cambia el inicio de traduccion
	$menor=$comienzo[0];
	for($i=1;$i<scalar(@comienzo);$i++){
		if($comienzo[$i]<$menor){
			$menor=$comienzo[$i];
		}
	}
####$menor tiene el valor mas bajo =P
####Tomar el valor de inicio de cada secuencia y concatenar espacios dependiendo del menor
	for($i=0;$i<$k;$i++){
		$dif=($comienzo[$i]-$menor)/3;
		$espacios='';
		for($j=0;$j<$dif-1;$j++){
			$espacios.=' ';
		}
		$proteina[$i]=$espacios.$proteina[$i];
	}

######
	$ref=$k-1;
	$cont=0;
	$penul="aas ref              ";
	$final="# aas identicos      ";
	for($i=0;$i<length($proteina[$ref]);$i++){
		$aaref=substr($proteina[$ref],$i,1);
		$simil=0;
		for($j=0;$j<$ref;$j++){
			$aacepa=substr($proteina[$j],$i,1);
			if($aacepa eq $aaref){
				$simil+=1;
			}
		}
		$penul.=$aaref."  ";
		$final.=$simil.' ';
		$finalito[$cont]=$simil;
		$cont+=1;		
	}
	#$output=">".$opts{'o'};
	$nombre=substr($file,0,length($file)-4);
	$output=">".$nombre."_prot.txt";
	$clustal=">".$nombre.".aln";
	open(OUT, $output);
	open(CLUSTAL, $clustal);
	print"Alineamiento del archivo $file traducido a aminoacidos\n\n";
	print OUT "Alineamiento del archivo $file traducido a aminoacidos\n\n";
	print CLUSTAL "CLUSTAL W 2.1 multiple sequence alignment\n\n";
	for($k=3;$k<3+$nseq;$k++){
		$cepas[$k-3]=substr($alin[$k],0,21);
		$cepas[$k-3].=$proteina[$k-3];
		print"$cepas[$k-3]\n";
		print OUT "$cepas[$k-3]\n";
		print CLUSTAL "$cepas[$k-3]\n";

	}
	$ncomp=$nseq-1;
	$asteriscos="                     ";        ####asumiendo que si se cambia el marco de lectura cambia toda la proteina
###############mejorando error
for($i=0;$i<scalar(@finalito);$i++){
	$naa2=$finalito[$i];
###	$anaa2=$finalito[$i-1];
###	#if($naa2<$anaa2){
	if($naa2<$nseq-1){
		$asteriscos.=" ";
	}
	else{
		$asteriscos.="*";
	}
}
############
#	for($i=22;$i<length($final);$i++){
#		$naa=substr($final,$i,1);
#		$anaa=substr($final,$i-1,1);
#		if($naa<$anaa){
#			$asteriscos.=" ";
#		}
#		else{
#			$asteriscos.="*";
#		}
#	}
	print"$asteriscos\n";
	print"$penul\n";
	print"$final\n";
	print OUT "$asteriscos\n";
	print CLUSTAL "$asteriscos\n";
	print OUT "$penul\n";
	print OUT "$final\n";
	close(OUT);
	close(CLUSTAL);
}
