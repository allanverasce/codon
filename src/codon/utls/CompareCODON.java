package codon.utls;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import codon.data.ORF;
import codon.recognizer.OrfRecognizer;

public class CompareCODON {
	
	public CompareCODON (String codon1, String codon2, String sameFile, String stopEquals1File, String stopEquals2File, String onlyIn1File, String onlyIn2File, String cds1File, String cds2File) {
		OrfRecognizer af1 = new OrfRecognizer ();
		OrfRecognizer af2 = new OrfRecognizer ();
		af1.loadFromCodonProject (new File (codon1));
		af2.loadFromCodonProject (new File (codon2));
		
		System.out.println (codon1 + " " + af1.candidateOrfCount ());
		System.out.println (codon2 + " " + af2.candidateOrfCount ());
		
		int sames = 0;
		int onlyIn1 = 0;
		int onlyIn2 = 0;
		int stopEquals = 0;
		String samesCoords = "";
		String stopEqualsStart1Coords = "";
		String stopEqualsStart2Coords = "";
		String onlyIn1Coords = "";
		String onlyIn2Coords = "";
		for (int i = 0; i < af1.orfsPorFrame.length; i++) {
			samesCoords += "FRAME " + i + "\n";
			onlyIn1Coords += "FRAME " + i + "\n";
			onlyIn2Coords += "FRAME " + i + "\n";
			for (ORF orf1: af1.orfsPorFrame [i]) {
				if (orf1.isRemoved () || orf1.product.entryId.equals ("RNA") || orf1.product.entryId.equals ("tRNA")) continue;
				boolean fnd = false;
				for (ORF orf2: af2.orfsPorFrame [i]) {
					if (orf2.isRemoved () || orf2.product.entryId.equals ("RNA") || orf2.product.entryId.equals ("tRNA")) continue;
					if (orf1.start == orf2.start && orf1.stop == orf2.stop) {
						sames ++;
						samesCoords += orf1.start + ".." + orf1.stop + ";" + orf1.product.gene + ";" + orf1.product.productName + ";" + orf1.getAccuracy () + "\n";
						fnd = true;
					}
					else if (orf1.isReverse () && (orf1.start == orf2.start) || (!orf1.isReverse ()) && (orf1.stop == orf2.stop)) {
						fnd = true;
						stopEqualsStart1Coords += orf1.start + ".." + orf1.stop + ";" + orf1.product.gene + ";" + orf1.product.productName + ";" + orf1.getAccuracy () + "\n";
						stopEqualsStart2Coords += orf2.start + ".." + orf2.stop + ";" + orf2.product.gene + ";" + orf2.product.productName + ";" + orf2.getAccuracy () + "\n";
						stopEquals ++;
					}
					//if (orf2.start > orf1.stop) break; 
				}
				if (!fnd) {
					onlyIn1 ++;
					onlyIn1Coords += orf1.start + ".." + orf1.stop + ";" + orf1.product.gene + ";" + orf1.product.productName + ";" + orf1.getAccuracy () + "\n";
				}
			}
			for (ORF orf2: af2.orfsPorFrame [i]) {
				if (orf2.isRemoved () || orf2.product.entryId.equals ("RNA") || orf2.product.entryId.equals ("tRNA")) continue;
				boolean fnd = false;
				for (ORF orf1: af1.orfsPorFrame [i]) {
					if (orf1.isRemoved () || orf1.product.entryId.equals ("RNA") || orf1.product.entryId.equals ("tRNA")) continue;
					if (orf1.start == orf2.start && orf1.stop == orf2.stop) {
						fnd = true;
					}
					else if (orf1.isReverse () && (orf1.start == orf2.start) || (!orf1.isReverse ()) && (orf1.stop == orf2.stop)) {
						fnd = true;
					}
					//if (orf1.start > orf2.stop) break;
				}
				if (!fnd) {
					onlyIn2 ++;
					onlyIn2Coords += orf2.start + ".." + orf2.stop + ";" + orf2.product.gene + ";" + orf2.product.productName + ";" + orf2.getAccuracy () + "\n";
				}
			}			
		}
		try {
			BufferedWriter writer = new BufferedWriter (new FileWriter (sameFile));
			writer.write (sames + " cds\n");
			writer.write (samesCoords);			
			writer.close ();
			writer = new BufferedWriter (new FileWriter (onlyIn1File));
			writer.write (onlyIn1 + " cds\n");
			writer.write (onlyIn1Coords);
			writer.close ();
			writer = new BufferedWriter (new FileWriter (onlyIn2File));
			writer.write (onlyIn2 + " cds\n");
			writer.write (onlyIn2Coords);
			writer.close ();
			writer = new BufferedWriter (new FileWriter (stopEquals1File));
			writer.write (stopEquals + " cds\n");
			writer.write (stopEqualsStart1Coords);
			writer.close ();
			writer = new BufferedWriter (new FileWriter (stopEquals2File));
			writer.write (stopEquals + " cds\n");
			writer.write (stopEqualsStart2Coords);
			writer.close ();
			writer = new BufferedWriter (new FileWriter (cds1File));
			int rnas = 0;
			int cds = 0;
			int genes = 0;
			int products = 0;
			int uncharacterized = 0;
			for (int i = 0; i < af1.orfsPorFrame.length; i++) {
				for (ORF orf1: af1.orfsPorFrame [i]) {
					if (orf1.product != null && orf1.product.entryId != null && (orf1.product.entryId.equals ("RNA") || orf1.product.entryId.equals ("tRNA"))) rnas ++;
					if (orf1.isRemoved () || orf1.product.entryId.equals ("RNA") || orf1.product.entryId.equals ("tRNA")) continue;
					if (orf1.isGene ()) genes ++;
					else if (orf1.isUncharacterized ()) uncharacterized ++;
					else products ++;
					cds ++;
					writer.write (">" + orf1.id);
					writer.newLine ();
					writer.write (orf1.getSequence ());
					writer.newLine ();
				}
			}
			writer.close ();
			System.out.println (cds1File + " " + rnas + " RNAS");
			System.out.println (cds1File + " " + genes + " GENES");
			System.out.println (cds1File + " " + products + " OTHER PRODUCT");
			System.out.println (cds1File + " " + uncharacterized + " UNCHARCTERIZED");
			System.out.println (cds1File + " " + cds + " cds");
			rnas = 0;
			cds = 0;
			genes = 0;
			products = 0;
			uncharacterized = 0;
			writer = new BufferedWriter (new FileWriter (cds2File));
			for (int i = 0; i < af2.orfsPorFrame.length; i++) {
				for (ORF orf2: af2.orfsPorFrame [i]) {
					if (!orf2.isRemoved () && orf2.product != null && orf2.product.entryId != null && (orf2.product.entryId.equals ("RNA") || orf2.product.entryId.equals ("tRNA"))) rnas ++;
					if (orf2.isRemoved () || orf2.product.entryId.equals ("RNA") || orf2.product.entryId.equals ("tRNA")) continue;
					if (orf2.isGene ()) genes ++;
					else if (orf2.isUncharacterized ()) uncharacterized ++;
					else products ++;
					cds ++;
					writer.write (">" + orf2.id);
					writer.newLine ();
					writer.write (orf2.getSequence ());
					writer.newLine ();
				}
			}
			writer.close ();
			System.out.println (cds2File + " " + rnas + " RNAS");
			System.out.println (cds2File + " " + genes + " GENES");
			System.out.println (cds2File + " " + products + " OTHER PRODUCT");
			System.out.println (cds2File + " " + uncharacterized + " UNCHARCTERIZED");
			System.out.println (cds2File + " " + cds + " cds");
			
		} catch (IOException e) {
			e.printStackTrace ();
		}
	}
	
	public static void main (String [] args) {
		new CompareCODON (
				"workspace/PALESB58/CODON.codon", 
				"workspace/PALESB58/NCBI_WITH_LOW_ACC.codon",
				"res/CN_EQUALS_CODON_NCBI.txt",
				"res/CN_SAME_STOP_CODON_products.txt",
				"res/CN_SAME_STOP_NCBI_products.txt",
				"res/CN_SPECIFICS_IN_CODON.txt",
				"res/CN_SPECIFICS_IN_NCBI.txt",
				"res/CDS_CODON.fasta",
				"res/CDS_NCBI.fasta");
		new CompareCODON (
				"workspace/PALESB58/CODON.codon", 
				"workspace/PALESB58/RAST_WITH_LOW_ACC.codon",
				"res/CR_EQUALS_CODON_RAST.txt",
				"res/CR_SAME_STOP_CODON_products.txt",
				"res/CR_SAME_STOP_RAST_products.txt",
				"res/CR_SPECIFICS_IN_CODON.txt",
				"res/CR_SPECIFICS_IN_RAST.txt",
				"res/CDS_CODON.fasta",
				"res/CDS_RAST.fasta");
		new CompareCODON (
				"workspace/PALESB58/NCBI_WITH_LOW_ACC.codon", 
				"workspace/PALESB58/RAST_WITH_LOW_ACC.codon",
				"res/RN_EQUALS_NCBI_RAST.txt",
				"res/RN_SAME_STOP_NCBI_products.txt",
				"res/RN_SAME_STOP_RAST_products.txt",
				"res/RN_SPECIFICS_IN_NCBI.txt",
				"res/RN_SPECIFICS_IN_RAST.txt",
				"res/CDS_NCBI.fasta",
				"res/CDS_RAST.fasta");
		new CompareCODON (
				"workspace/PALESB58/CODON.codon", 
				"workspace/PALESB58/NCBI_WITHOUT_LOW_ACC.codon",
				"res/CNO_EQUALS_CODON_NCBI.txt",
				"res/CNO_SAME_STOP_CODON_products.txt",
				"res/CNO_SAME_STOP_NCBI_products.txt",
				"res/CNO_SPECIFICS_IN_CODON.txt",
				"res/CNO_SPECIFICS_IN_NCBI.txt",
				"res/CDS_CODON.fasta",
				"res/CDS_NCBI_O.fasta");
		new CompareCODON (
				"workspace/PALESB58/CODON.codon", 
				"workspace/PALESB58/RAST_WITHOUT_LOW_ACC.codon",
				"res/CRO_EQUALS_CODON_RAST.txt",
				"res/CRO_SAME_STOP_CODON_products.txt",
				"res/CRO_SAME_STOP_RAST_products.txt",
				"res/CRO_SPECIFICS_IN_CODON.txt",
				"res/CRO_SPECIFICS_IN_RAST.txt",
				"res/CDS_CODON.fasta",
				"res/CDS_RAST_O.fasta");
		new CompareCODON (
				"workspace/PALESB58/NCBI_WITHOUT_LOW_ACC.codon", 
				"workspace/PALESB58/RAST_WITHOUT_LOW_ACC.codon",
				"res/RNO_EQUALS_NCBI_RAST.txt",
				"res/RNO_SAME_STOP_NCBI_products.txt",
				"res/RNO_SAME_STOP_RAST_products.txt",
				"res/RNO_SPECIFICS_IN_NCBI.txt",
				"res/RNO_SPECIFICS_IN_RAST.txt",
				"res/CDS_NCBI_O.fasta",
				"res/CDS_RAST_O.fasta");
	}

}
