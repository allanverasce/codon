package codon.data;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Vector;

public class OrganismFamilly extends ArrayList <Organism> {
	private static final long serialVersionUID = 1L;
	
	public String famillyName;
	public int occurrences = 0;
	public int genes = 0;
	public int proteins = 0;
	public int uncharacterized = 0;
	public Vector <ORF> orfs = new Vector <ORF> ();
	
	public OrganismFamilly (String organism) {
		this.famillyName = getFamillyName (organism);
	}
	
	public void sort () {
		sort (new Comparator <Organism> () {
			public int compare (Organism org0, Organism org1) {
				if (org0.occurrences > org1.occurrences) return -1;
				if (org0.occurrences < org1.occurrences) return 1;
				return 0;
			}
		});
	}
	
	public boolean isFromSameFamilly (String organism) {
		return famillyName.equals (getFamillyName (organism));
	}
	
	public static String getFamillyName (String organism) {
		String familly = organism.contains (" ") ? organism.substring (0, organism.indexOf (" ")) : organism;
		return familly;
	}
	
	public void add (ORF orf) {
		orfs.add (orf);
		occurrences ++;
		if (orf.isGene ()) genes ++;
		else if (orf.isUncharacterized ()) uncharacterized ++;
		else proteins ++;
		for (int i = 0; i < size (); i++) {
			if (get (i).organismName.equals (orf.product.organism)) {
				get (i).addOrf (orf);
				return;
			}
		}
		Organism org = new Organism (orf.product.organism);
		org.addOrf (orf);
		add (org);
	}
}
