package codon.utls;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class BlastInfoAnalyzer {

	class BlastInfoAnalysisResult {
		String id = "";
		String entry = "";
		boolean removed = false;
		int initialLength = 0;
		int finalLength = 0;
		boolean blasted = false;
		int killed = 0;
		boolean gene = false;
		boolean uncharacherized = false;
		float accuracy = 0;
		public String toString () {
			return id + ";" + entry + ";" + removed + ";" + initialLength + ";" + finalLength + ";" + blasted + ";" + killed + ";" + gene + ";" + uncharacherized + ";" + accuracy;
		}
	}
	
	BlastInfoAnalysisResult loadBlastInfoAnalysisResult (String s) {
		BlastInfoAnalysisResult res = new BlastInfoAnalysisResult ();
		String items [] = s.split (";");
		res.id = items [0];
		res.entry = items [1];
		res.removed = Boolean.parseBoolean (items [2]);
		res.initialLength = Integer.parseInt (items [3]);
		res.finalLength = Integer.parseInt (items [4]);
		res.blasted = Boolean.parseBoolean (items [5]);
		res.killed = Integer.parseInt (items [6]);
		res.gene = Boolean.parseBoolean (items [7]);
		res.uncharacherized = Boolean.parseBoolean (items [8]);
		res.accuracy = Float.parseFloat (items [9].replace (",", "."));
		return res;
	}
	
	ArrayList <BlastInfoAnalysisResult> results = new ArrayList <BlastInfoAnalyzer.BlastInfoAnalysisResult> ();
	
	public BlastInfoAnalyzer (ArrayList <String> files) {
		for (String file: files) {
			try {
				BufferedReader reader = new BufferedReader (new FileReader (file));
				reader.readLine ();
				String line;
				while ((line = reader.readLine ()) != null) {
					results.add (loadBlastInfoAnalysisResult (line));
				}
				reader.close ();
			} catch (IOException e) {
				e.printStackTrace ();
			}
		}
	}
	
	public boolean isIn (BlastInfoAnalysisResult res, ArrayList <String> filters) {
		boolean isIn = true;
		for (String f : filters) {
			if (f.indexOf ("I") == 0) {
				int val = Integer.parseInt (f.substring (1));
				if (res.initialLength < val) isIn = false; 
			}
			else if (f.indexOf ("i") == 0) {
				int val = Integer.parseInt (f.substring (1));
				if (res.initialLength > val) isIn = false; 
			}
			else if (f.indexOf ("F") == 0) {
				int val = Integer.parseInt (f.substring (1));
				if (res.finalLength < val) isIn = false; 
			}
			else if (f.indexOf ("f") == 0) {
				int val = Integer.parseInt (f.substring (1));
				if (res.finalLength > val) isIn = false; 
			}
			else if (f.indexOf ("K") == 0) {
				int val = Integer.parseInt (f.substring (1));
				if (res.killed < val) isIn = false; 
			}
			else if (f.indexOf ("k") == 0) {
				int val = Integer.parseInt (f.substring (1));
				if (res.killed > val) isIn = false; 
			}
			else if (f.indexOf ("A") == 0) {
				float val = Float.parseFloat (f.substring (1).replace (",", "."));
				if (res.accuracy < val) isIn = false; 
			}
			else if (f.indexOf ("a") == 0) {
				float val = Float.parseFloat (f.substring (1).replace (",", "."));
				if (res.accuracy > val) isIn = false;
			}
			else if (f.indexOf ("B") == 0) {
				if (!res.blasted) isIn = false; 
			}
			else if (f.indexOf ("b") == 0) {
				if (res.blasted) isIn = false; 
			}
			else if (f.indexOf ("B") == 0) {
				if (!res.blasted) isIn = false; 
			}
			else if (f.indexOf ("b") == 0) {
				if (res.blasted) isIn = false; 
			}
			else if (f.indexOf ("G") == 0) {
				if (!res.gene) isIn = false; 
			}
			else if (f.indexOf ("g") == 0) {
				if (res.gene) isIn = false; 
			}
			else if (f.indexOf ("D") == 0) {
				if (!res.removed) isIn = false; 
			}
			else if (f.indexOf ("d") == 0) {
				if (res.removed) isIn = false; 
			}
			else if (f.indexOf ("U") == 0) {
				if (!res.uncharacherized) isIn = false; 
			}
			else if (f.indexOf ("u") == 0) {
				if (res.uncharacherized) isIn = false; 
			}
			else if (f.indexOf ("R") == 0) {
				if (!res.entry.equals ("RNA") && !res.entry.equals ("tRNA")) isIn = false; 
			}
			else if (f.indexOf ("r") == 0) {
				if (res.entry.equals ("RNA") || res.entry.equals ("tRNA")) isIn = false; 
			}
		}
		return isIn;
	}
	
	public int countFilter (ArrayList <String> filters) {
		int count = 0;
		for (BlastInfoAnalysisResult res: results) {
			if (isIn (res, filters)) count ++;
		}
		return count;
	}
	
	public void list (ArrayList <String> filters) {
		int count = 0;
		float accuracy = 0;
		float initialL = 0;
		float finalL = 0;
		int genes = 0;
		int products = 0;
		int uncharacherized = 0;
		for (BlastInfoAnalysisResult res: results) {
			if (isIn (res, filters)) {
				count ++;
				System.out.println (res);
				accuracy += res.accuracy;
				initialL += res.initialLength;
				finalL += res.finalLength;
				if (res.gene) genes ++;
				else if (!res.uncharacherized) products ++;
				else uncharacherized ++;
			}
		}
		accuracy /= count;
		initialL /= count;
		finalL /= count;
		System.out.println (count + " results, genes:" + genes + " other-protein:" + products + " uncharacterized-protein:" + uncharacherized + "\naverages: accuracy=" + accuracy  + "  initial-lengh=" + initialL + "  final-length" + finalL);
	}
	
	public ArrayList <String> createFilter (String filter) {
		String fi [] = filter.split (";");
		ArrayList <String> fs = new ArrayList <> ();
		for (String s: fi) fs.add (s);
		return fs;
	}
	
	public static void main (String [] args) {
		ArrayList <String> files = new ArrayList <> ();
//		files.add ("C:\\workspace\\NCTC11134.csv");
//		files.add ("C:\\workspace\\CP026353.1.csv");
		files.add ("C:\\workspace\\LESB58.csv");
//		files.add ("C:\\workspace\\KP5-1.csv");
//		files.add ("C:\\workspace\\H37Rv.csv");
//		files.add ("C:\\workspace\\CP041647.1.csv");
//		files.add ("C:\\workspace\\M26365.csv");
//		files.add ("C:\\workspace\\FDAARGOS_768.csv");
//		files.add ("C:\\workspace\\NZ_AP019695-FA.1.csv");
		BlastInfoAnalyzer analizer = new BlastInfoAnalyzer (files);
		ArrayList <String> filter1 = analizer.createFilter ("i90;f90;A80;d;r;u");
		analizer.list (filter1);
//		filter1 = analizer.createFilter ("B");
//		analizer.list (filter1);
//		filter1 = analizer.createFilter ("i110;f110;A80;u;d;G");
//		analizer.list (filter1);
//		filter1 = analizer.createFilter ("i120;f120;A80;u;d;G");
//		analizer.list (filter1);
	}
}
