package codon.blast;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.Vector;

import javax.swing.Timer;

import codon.data.CodonConfiguration;
import codon.data.ORF;
import codon.recognizer.OrfRecognizer;
import codon.ui.PlotOrfsFrame;
import splash.WaitingSplash;
import splash.WaitingSplashThread;

public class BlastingManager implements BlastingThreadListener {

	Vector <BlastingThread> threads = new Vector <> ();
	int threadNb = 1;
	
	OrfRecognizer reco;
	List <ORF> validOrfs = new ArrayList 	<> ();
	
	public int previouslyDone = 0;
	public int ignored = 0;
	public int done = 0;
	boolean fullBlast = false;
	public int stopped = 0;
	
	ArrayList <BlastingManagerListener> listeners = new ArrayList <BlastingManagerListener> ();
	
	public void addBlastingManagerListener (BlastingManagerListener l) {
		listeners.add (l);
	}
	
	public void removeBlastingManagerListener (BlastingManagerListener l) {
		listeners.remove (l);
	}
	
	public void fireThreadStatustChanged (BlastingThread blastingThread) {
		for (BlastingManagerListener l: listeners) l.threadStatustChanged (this, blastingThread); 
	}
	
	public void fireEnded (BlastingThread endedThread) {
		for (BlastingManagerListener l: listeners) l.threadEnded (this, endedThread); 
	}
	
	public void firePreviouslyDone (ORF orf, int toProcess) {
		for (BlastingManagerListener l: listeners) l.orfPreviouslyDone (this, orf, toProcess); 
	}
	
	public void fireIgnored (ORF orf) {
		for (BlastingManagerListener l: listeners) l.orfIgnored (this, orf); 
	}
	
	public void fireReincluded (ORF orf) {
		for (BlastingManagerListener l: listeners) l.orfReincluded (this, orf); 
	}
	
	public void fireRebounded (ORF orf) {
		for (BlastingManagerListener l: listeners) l.orfRebounded (this, orf, orf.start, orf.stop); 
	}
	
	public void fireReboundAndRename (ORF orf, int newStart, int newStop) {
		for (BlastingManagerListener l: listeners) l.orfReboundedAndRename (this, orf, newStart, newStop); 
	}
	
	public void fireStarted (BlastingThread endedThread) {
		for (BlastingManagerListener l: listeners) l.threadStarted (this, endedThread); 
	}
	
	public void fireRestarted (BlastingThread endedThread) {
		for (BlastingManagerListener l: listeners) l.threadRestarted (this, endedThread); 
	}
	
	public void fireFinalized () {
		for (BlastingManagerListener l: listeners) l.blastingFinalized (this); 
	}
	
	class ThreadTimerAction implements ActionListener {
		BlastingThread blastingThread;
		public ThreadTimerAction (BlastingThread blastingThread) {
			this.blastingThread = blastingThread;
		}
		@SuppressWarnings ("deprecation")
		public void actionPerformed (ActionEvent e) {
			Timer tm = (Timer) e.getSource ();
			if (blastingThread.blaster.receiving) {
				blastingThread.blaster.receiving = false;
				tm.restart ();
				return;
			}
			blastingThread.stopped = true;
			stopped ++;
			blastingThread.timer.stop ();
//			blastingThread.interrupt ();
			blastingThread.stop ();
			threads.remove (blastingThread);
			BlastingThread th  = new BlastingThread (blastingThread.orf, reco.dataDestDirS);
			fireRestarted (th);
			th.addBlastingThreadListener (BlastingManager.this);
			tm = new Timer (getTimeOut (th.orf), new ThreadTimerAction (th));
			th.timer = tm;
			th.start ();
			threads.add (th);
			tm.start ();
		}
	}
	
	LinkedList <Integer> lastTimes = new LinkedList <> ();
	
	public void insertTime (int t) {
		lastTimes.add (t);
		if (lastTimes.size () > 20) lastTimes.remove (0);
	}
	
	public int getTimeOut (ORF orf) {
		if (lastTimes.size () < 20) return orf.length * 4 * 1000;  
		else {
			int avg = 0;
			for (int i = 0; i < lastTimes.size (); i++) {
				avg += lastTimes.get (i);
			}
			avg /= lastTimes.size ();
			return 3 * avg;
		}
	}
	
	public BlastingManager (OrfRecognizer reco, int threadNumber, boolean fullBlast) {
		this.reco = reco;
		this.fullBlast = fullBlast;
		this.threadNb = threadNumber;
	}

	WaitingSplash splash;
	public void start (List <ORF> validOrfs) {
		for (List <ORF> orfs: reco.orfsPorFrame) {
			for (ORF orf: orfs) {
				orf.marked = orf.isRemoved ();
			}
		}
		this.validOrfs = validOrfs;
		validOrfs.sort (new Comparator <ORF> () {
			public int compare (ORF orf1, ORF orf2) {
				if (orf1.length > orf2.length) return 1;
				if (orf1.length == orf2.length) return 0;
				return -1;
			}
		});
		WaitingSplashThread th = new WaitingSplashThread () {
			public void run () {
				try {
					Thread.sleep (100);
				} catch (InterruptedException e) {
					e.printStackTrace 	();
				}
				for (int j = 0; j < threadNb; j++) {
					doBlast ();
				}
				if (splash != null) splash.setVisible (false);
			}
		};
		if (CodonConfiguration.ui) {
			try {
				Thread.sleep (100);
			} catch (InterruptedException e) {
				e.printStackTrace 	();
			}
			splash = new WaitingSplash (PlotOrfsFrame.instance, "Processing previous blasted ORF", "No predicted time to end... It will depend on how many ORFs have ever been blasted. Wait! ^_^", 100, th, false);
			th.start ();
			splash.setVisible (true);
		}
		else {
			th.run ();
		}
		if (toProcess () == 0 && threads.size () == 0) fireFinalized ();
	}
	
	@SuppressWarnings("deprecation")
	public void interrupt () {
		for (int i = 0; i < threads.size (); i++) {
			threads.get (i).timer.stop ();
			threads.get (i).stopped = true;
			threads.get (i).stop ();
		}
		threads.clear ();
	}
	
	synchronized public boolean doBlast () {
		ORF orf;
		int prev = 0;
		while (validOrfs.size () != 0) {
			orf = validOrfs.get (validOrfs.size () - 1);
			orf.reloadProduct ();
			if (orf.isBlasted () || orf.hasSwissUpdate ()) {
				validOrfs.remove (validOrfs.size () - 1);
				if (!orf.marked) {
					prev ++;
					previouslyDone ++;
					firePreviouslyDone (orf, validOrfs.size ());
				}
				orf.marked = true;
				if (!fullBlast) {
					orf.autoRebound ();
					if (validOrfs.size () == 0 || orf.length >= validOrfs.get (validOrfs.size () - 1).length) cleanOverlap (orf);
					else {
						reinsertOrf (orf);
						fireRebounded (orf);
					}
				}
				if (prev % 100 == 0 && splash != null && splash.isVisible ()) splash.changeStatus ("No predicted time to end... It will depend on how many ORFs have ever been blasted. Wait! ^_^ (did: " + prev + ")", (int) (Math.random () * 100));
			}
			else {
				break;
			}
		}
		if (validOrfs.size () == 0) return false;
		else orf = validOrfs.remove (validOrfs.size () - 1);
		orf.marked = true;
		BlastingThread th = new BlastingThread (orf, reco.dataDestDirS);
		th.addBlastingThreadListener (this);
		fireStarted (th);
		threads.add (th);
		Timer tm = new Timer (getTimeOut (th.orf), new ThreadTimerAction (th));
		th.timer = tm;
		tm.start ();
		th.start ();
		return true;
	}
	
	public int toProcess () {
		return validOrfs.size ();
	}
	
	public int ongoing () {
		return threads.size ();
	}
	
	private void cleanOverlap (ORF mainOrf) {
		mainOrf.reloadProduct ();
		mainOrf.autoRebound ();
//		mainOrf.selectBetterEntryMaximizingIntergenicRegion ();
		if (mainOrf.id.equals ("F2-000936-2832340-2833276")) {
			System.err.println (mainOrf);
			System.err.println (mainOrf.getAccuracy () + "  " + mainOrf.sizesFit () + "   " + mainOrf.isGene ());
		}
		boolean redo = false;
		if (mainOrf.getAccuracy () > CodonConfiguration.accuracy_level_to_remove_overlap_during_blast && mainOrf.sizesFit ()) {
			for (List <ORF> orfs: reco.orfsPorFrame) {
				for (ORF orf: orfs) {
					if (orf.isRemoved () || !mainOrf.isGene () && orf.length > mainOrf.length / 2) continue;
					int intersection = orf.intersection (mainOrf);
					if (mainOrf != orf && intersection > 0 
							&&
							(!orf.isReverse () && !mainOrf.isReverse () && orf.stop > mainOrf.stop
							|| orf.isReverse () && mainOrf.isReverse () && orf.start < mainOrf.start
							|| orf.isReverse () && !mainOrf.isReverse () && orf.start < mainOrf.start
							|| !orf.isReverse () && mainOrf.isReverse () && orf.stop > mainOrf.stop) 
							&& 
							!orf.marked) {
							int queue = orf.intersection (mainOrf);
						if (!orf.isReverse () && !mainOrf.isReverse () && orf.start < mainOrf.start) queue += mainOrf.start - orf.start;
						else if (orf.isReverse () && mainOrf.isReverse () && orf.stop > mainOrf.stop) queue += orf.stop - mainOrf.stop;
						else if (orf.isReverse () && !mainOrf.isReverse () && orf.stop > mainOrf.stop) queue += orf.stop - mainOrf.stop;
						else if (!orf.isReverse () && mainOrf.isReverse () && orf.start < mainOrf.start) queue += mainOrf.start - orf.start;
						orf.castrating (queue);
						if (orf.length < CodonConfiguration.min_orf_length) {
							orf.setRemoved (true);
						}
						else {
							fireReboundAndRename (orf, orf.start, orf.stop);
							orf.setOriginalStart (orf.start);
							orf.setOriginalStop (orf.stop);
							orf.rename ();
							orf.reloadProduct ();
							validOrfs.remove (orf);
							if (!orf.isRemoved ()) reinsertOrf (orf);
						}
					}
					else if (mainOrf != orf && intersection > CodonConfiguration.overlap_auto_removed * orf.length / 100) {
						if (orf.marked) {
							if (orf.getAccuracy () < 80) orf.setRemoved (true);
							else continue;
						}
						else {
							orf.setRemoved (true);
						}
					}
					if (orf.isRemoved ()) {
						mainOrf.killed.add (orf);
						for (ORF killed: orf.killed) {
							killed.marked = false;
							ignored --;
							fireReincluded (killed);
							reinsertOrf (killed);
							killed.setRemoved (false);
							redo = true;
						}
						orf.marked = true;
						validOrfs.remove (orf);
						ignored ++;
						fireIgnored (orf);
					}
				}
			}
		}
		if (redo) cleanOverlap (mainOrf);
	}
	
	private void reinsertOrf (ORF orf) {
		boolean fnd = false;
		for (int i = 0; i < validOrfs.size (); i ++) {
			if (validOrfs.get (i).length > orf.length) {
				fnd = true;
				validOrfs.add (i, orf);
				break;
			}
		}
		if (!fnd) validOrfs.add (orf);
	}

	synchronized public void addOrfToBlast (ORF orf) {
		validOrfs.add (orf);
		if (threads.size () < threadNb) doBlast ();
	}
	
	public void threadStatustChanged (BlastingThread blastingThread) {
		fireThreadStatustChanged (blastingThread);
	}

	public void threadEnded (BlastingThread blastingThread) {
		ended (blastingThread);
	}
	
	synchronized public void ended (BlastingThread thread) {
		thread.timer.stop ();
		threads.remove (thread);
		done ++;
		insertTime ((int) (System.currentTimeMillis () - thread.startTime));
		ORF mainOrf = thread.orf;
		if (!fullBlast) {
			if (validOrfs.size () == 0 || mainOrf.length >= validOrfs.get (validOrfs.size () - 1).length) cleanOverlap (mainOrf);
			else reinsertOrf (mainOrf);
		}
		if (threads.size () < (int) threadNb) doBlast ();
		fireEnded (thread);
		if (toProcess () == 0 && threads.size () == 0) fireFinalized ();
	}

	public void setMaxThread (int nb) {
		for (int i = threads.size (); i < nb; i ++) doBlast ();
		threadNb = nb;
	}
		
	@SuppressWarnings("deprecation")
	synchronized public void kill (int ind) {
		BlastingThread thread = threads.get (ind);
		thread.stopped = true;
		stopped ++;
		thread.timer.stop ();
		thread.stop ();
		threads.remove (thread);
		reinsertOrf (thread.orf);
		setMaxThread (threadNb); 
		//if (threads.size () < (int) threadNb) doBlast ();
	}

	public void restart (String orfId) {
		for (int i = 0; i < threads.size (); i++) {
			BlastingThread thread = threads.get (i); 
			if (thread.orf.id.equals (orfId)) {
				kill (i);
				break;
			}
		}
	}
	
}
