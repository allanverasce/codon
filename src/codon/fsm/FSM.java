package codon.fsm;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Hashtable;

public class FSM implements Serializable {
	private static final long serialVersionUID = 1L;
	static char [] bases = new char [] {'A', 'T', 'G', 'C'};
	
	class Transition implements Serializable {
		private static final long serialVersionUID = 1L;
		char read;
		FSMState nextState;
		public Transition (char read, FSMState nextState) {
			super ();
			this.read = read;
			this.nextState = nextState;
		}
	}
	
	class FSMState implements Serializable {
		private static final long serialVersionUID = 1L;
		ArrayList <Transition> transistions = new ArrayList <> ();
		String recognized = null;
		int length = -1;
		FSMState addTransition (char c, String prefix) {
			for (Transition tr: transistions) {
				if (tr.read == c) return tr.nextState;
			}
			Transition tr = new Transition (c, new FSMState ());
			transistions.add (tr);
			return tr.nextState;
		}
		FSMState next (char c) {
			for (Transition tr: transistions) {
				if (tr.read == c) return tr.nextState;
			}
			return null;
		}
	}
	
	FSMState idleState = new FSMState ();
	
	public void mark (String src) {
		FSMState current = idleState; 
		for (int i = 0; i < src.length (); i++) {
			char c = src.charAt (i);
			current = current.next (c);
			if (current.recognized != null) System.err.println ("> " + current.recognized);
		}
	}
	
	public boolean load (String filename) {
		Hashtable <String, FSMState> prefixs = new Hashtable <String, FSM.FSMState> ();
		try {
			File file = new File (filename);
			if (!file.exists ()) return false;
			BufferedReader reader = new BufferedReader (new FileReader (file));
			int nb = Integer.parseInt (reader.readLine ());
			for (int i = 0; i < nb; i++) {
				String pref = reader.readLine ();
				prefixs.put (pref, new FSMState ());
			}
			idleState = prefixs.get ("");
			for (int i = 0; i < nb; i++) {
				String l = reader.readLine ();
				String trans [] =  l.split ("\\|");
				
				for (int j = 0; j < trans.length; j += 4) {
					FSMState current = prefixs.get (trans [j]);
					FSMState next = prefixs.get (j + 2 >= trans.length ? "" : trans [j + 2]);
					Transition tr = new Transition (trans [j + 1].charAt (0), next);
					if (j + 3 < trans.length && !trans [j + 3].equals ("")) {
						next.recognized = trans [j + 3];
						next.length = Integer.parseInt (trans [j + 4]);
					}
					current.transistions.add (tr);
				}
			}
			reader.close ();
		}
		catch (IOException e) {
			e.printStackTrace ();
		} 
		return true;
	}
	
}
