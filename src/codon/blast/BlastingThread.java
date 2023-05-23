package codon.blast;

import java.util.ArrayList;

import javax.swing.Timer;

import codon.data.ORF;

public class BlastingThread extends Thread {
	Timer timer;
	public boolean withoutConnection = false;
	boolean stopped = false;
	
	Blaster blaster = new Blaster ();
	long startTime = 0; 
	
	public ORF orf = null;
	String dataDestFileS = null;
	
	public BlastingThread (ORF toPtocess, String dataDestFileS) {
		orf = toPtocess;
		this.dataDestFileS = dataDestFileS;
	}
	
	ArrayList <BlastingThreadListener> listeners = new ArrayList <> ();
	public void addBlastingThreadListener (BlastingThreadListener l) {
		listeners.add (l);
	}
	public void removeBlastingThreadListener (BlastingThreadListener l) {
		listeners.remove (l);
	}
	public void fireThreadStatustChanged () {
		for (BlastingThreadListener l: listeners) l.threadStatustChanged (this); 
	}
	public void fireEnded () {
		for (BlastingThreadListener l: listeners) l.threadEnded (this); 
	}
	
	public void run () {
		startTime = System.currentTimeMillis ();
		while (true) {
			if (blaster.blast (orf, dataDestFileS) != -1) {
				orf.setLoaded (false);
				orf.reloadProduct ();
				orf.autoRebound ();
				if (stopped) return;
				timer.restart ();
				break;
			}
			else {
				withoutConnection = true;
				fireThreadStatustChanged ();
				while (!InternetConnectivityTester.isConnected ()) {
					try {
						if (stopped) return;
						timer.restart ();
						Thread.sleep (30000);
						
					} catch (InterruptedException e1) {
						e1.printStackTrace ();
					}
				}
				withoutConnection = false;
				fireThreadStatustChanged ();
			}
		}
		fireEnded ();
	}
	
	public String toString () {
		String sec = (System.currentTimeMillis () - startTime) / 1000 + "";
		while (sec.length () < 4) sec = "0" + sec;
		return  sec + "s " + (withoutConnection ? "(paused) " : " ") + (blaster.swissOngoing ? "S ": "U ")  + orf.id;
	}
}
