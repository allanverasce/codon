package codon.data;

public class GeneFeature {

	public String gene;
	public int start;
	public int stop;
	public boolean reverse;
	public boolean pseudo;
	public GeneFeature (String gene, int start, int stop, boolean reverse, boolean pseudo) {
		this.gene = gene;
		this.start = start;
		this.stop = stop;
		this.reverse = reverse;
		this.pseudo = pseudo;
	}
	
	public GeneFeature (int start, int stop, boolean reverse, boolean pseudo) {
		this.start = start;
		this.stop = stop;
		this.reverse = reverse;
		this.pseudo = pseudo;
	}

	public int getFrame () {
		return (start - 2 ) % 3 + (reverse ? 3: 0);
	}
}
