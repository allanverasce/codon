package codon.data;

import java.util.ArrayList;
import java.util.Hashtable;

public class LexicographicTree {

	ArrayList <String> duplicated = new ArrayList <> ();
	
	class LexicoNode {
		char root;
		Hashtable <Character, LexicoNode> childs = new Hashtable <Character, LexicoNode> ();
		ArrayList <ORF> orfs = new ArrayList <> ();
		public LexicoNode (char root) {
			this.root = root;
		}
		public void addString (String ref, String s, ORF orf) {
			if (s.length() == 0) {
				orfs.add (orf);
				if (orfs.size () == 2) duplicated.add (ref);
				return;
			}
			char c = s.charAt (0);
			LexicoNode nd = childs.get (c);
			if (nd == null) {
				nd = new LexicoNode (c);
				childs.put (c, nd);
			}
			nd.addString (ref, s.substring (1), orf);
		}
		public ArrayList <ORF> get (String s) {
			if (s.length () == 0) return orfs;
			char c = s.charAt (0);
			LexicoNode nd = childs.get (c);
			if (nd == null) return null;
			return nd.get (s.substring (1));
		}
		public boolean isDuplicated (String s) {
			if (s.length () == 0) return orfs.size () > 1;
			char c = s.charAt (0);
			LexicoNode nd = childs.get (c);
			if (nd == null) return false;
			return nd.isDuplicated (s.substring (1));
		}
	}
	
	LexicoNode root = new LexicoNode ('*');
	
	public void addString (String ref, ORF orf) {
		root.addString (ref, ref, orf);
	}
	
	public ArrayList <ORF> get (String s) {
		return root.get (s);
	}
	
	public boolean isDuplicated (String s) {
		return root.isDuplicated (s);
	}
}
