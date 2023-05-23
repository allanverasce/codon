package codon.data;

public class Codon {
	public int type;
	public int index;
	public Codon (int type, int index) {
		this.type = type;
		this.index = index;
	}
	public String toString () {
		return type + "º" + index;
	}
	public static Codon createFromString (String line) {
		String arg [] = line.split ("º");
		return new Codon (Integer.parseInt (arg [0]), Integer.parseInt (arg [1]));
	}
}
