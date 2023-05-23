package codon.blast;

public interface BlastingMonitorListener {
	public void blastRestarted (BlastingMonitor manager);
	public void blastRestarted (BlastingMonitor manager, String orfId);
	public void blastStarted (BlastingMonitor manager, String orfId);
	public void blastEnded (BlastingMonitor manager, String orfId);
	public void orfIgnored (BlastingMonitor manager, String orfId);
	public void orfRebounded (BlastingMonitor manager, String orfId, int start, int stop);
	public void orfReboundedAndRenamed (BlastingMonitor manager, String orfId, int start, int stop);
	public void blastPreviouslyDone (BlastingMonitor manager, String orfId);
	public void blastingFinalized (BlastingMonitor manager);
	public void orfReincluded (BlastingMonitor manager, String orf);
	public void statusChanged (BlastingMonitor blastingMonitor, String orfId, boolean status);
}
