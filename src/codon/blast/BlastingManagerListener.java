package codon.blast;

import codon.data.ORF;

public interface BlastingManagerListener {
	public void threadStatustChanged (BlastingManager manager, BlastingThread endedThread);
	public void threadStarted (BlastingManager manager, BlastingThread endedThread);
	public void threadRestarted (BlastingManager manager, BlastingThread endedThread);
	public void threadEnded (BlastingManager manager, BlastingThread endedThread);
	public void orfIgnored (BlastingManager manager, ORF orf);
	public void orfReincluded (BlastingManager manager, ORF orf);
	public void orfRebounded (BlastingManager manager, ORF orf, int start, int stop);
	public void orfReboundedAndRename (BlastingManager manager, ORF orf, int start, int stop);
	public void orfPreviouslyDone (BlastingManager manager, ORF orf, int toProcess);
	public void blastingFinalized (BlastingManager manager);
}
