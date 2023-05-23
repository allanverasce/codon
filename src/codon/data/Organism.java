package codon.data;

import java.util.Vector;

public class Organism {
	public String organismName = "";
	public int occurrences = 0;
	public int genes = 0;
	public int proteins = 0;
	public int uncharacterized = 0;
	
	public Vector <ORF> orfs = new Vector <> ();
	
	public Organism (String organism) {
		this.organismName = organism;
	}
	
	public void addOrf (ORF orf) {
		occurrences ++;
		if (orf.isGene ()) genes ++;
		else if (orf.isUncharacterized ()) uncharacterized ++;
		else proteins ++;
		orfs.add (orf);
	}
}
