package codon.fsm;

public class RecognizedPattern {
	
	public String recognized;
	public int start;
	public int end;

	public RecognizedPattern (String recognized, int start, int end) {
		this.recognized = recognized;
		this.start = start;
		this.end = end;
	}

	public int length () {
		return end - start;
	}
	
	public String toString () {
		return start + ";" + end + ";" + recognized ;
	}
}
