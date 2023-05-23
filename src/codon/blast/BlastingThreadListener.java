package codon.blast;

public interface BlastingThreadListener {
	public void threadStatustChanged (BlastingThread blastingThread);
	public void threadEnded (BlastingThread blastingThread);
}
