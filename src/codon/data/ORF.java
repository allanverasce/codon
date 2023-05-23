package codon.data;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.text.similarity.LevenshteinDetailedDistance;
import org.apache.commons.text.similarity.LevenshteinResults;

import codon.recognizer.OrfRecognizer;

public class ORF {
	
//	public final static int MIN_LENGHT = 80;
	
	public int length;
	public int start;
	public int stop;
	public int startOriginal;
	public int stopOriginal;
	private boolean removed = false;
	
	private boolean excluded = false;
	private boolean forced = false;
	public String id;
	public ProductProtein product;
	public OrfRecognizer reco = null;
	public int frame;
	
	public boolean marked = false;
	
	public static LevenshteinDetailedDistance levDDist = new LevenshteinDetailedDistance ();
	
	public ArrayList <ORF> killed = new ArrayList <> ();
	
	public ORF (OrfRecognizer reco, int start, int stop, boolean reverse) {
		this (reco, start, stop, (start + 1) % 3 + 3 * (reverse ? 1 : 0));
	}
	
	public ORF (OrfRecognizer reco, int start, int stop, int frame) {
		this.length = Math.abs (start - stop);
		if (length < CodonConfiguration.min_orf_length) setRemoved (true);
		this.start = start;
		this.stop = stop;
		this.startOriginal = start;
		this.stopOriginal = stop;
		this.frame = frame;
		this.reco = reco;
		rename ();
	}
	
	public void rename () {
		String lg = "" + length;
		for (int k = lg.length (); k < 6; k ++) lg = '0' + lg; 
		id = "F" + frame + "-" + lg + "-" + startOriginal + "-" + stopOriginal;
	}
	
	public void setOriginalStart (int start) {
		this.start = start;
		this.startOriginal = start;
		this.length = stop - start;
		setRemoved (length < CodonConfiguration.min_orf_length);
	}
	
	public void setOriginalStop (int stop) {
		this.stop = stop;
		this.stopOriginal = stop;
		this.length = stop - start;
		setRemoved (length < CodonConfiguration.min_orf_length);
	}
	
	public int intersection (ORF p) {
		return intersection (p.start, p.stop);
	}
	
	public int intersection (int begin, int end) {
		if (begin > stop || start > end) return 0;
		if (start <= begin && begin < stop) {
			if (end <= stop) return (end - begin);
			else return stop - begin;
		}
		if (begin < start && start < end) {
			if (stop <= end) return length;
			else return end - start;
		}
		return 0;
	}
	
	public int [] getQueryMatchBounds () {
		return getQueryMatchBounds (product);
	}
	
	public int [] getQueryMatchBounds (ProductProtein p) {
		int oldStartQuery = p.getStartQuerySeq ();
		int oldEndQuery = p.getEndQuerySeq ();
		int oldStartMatch = p.getStartMatchSeq ();
		int oldEndMatch = p.getEndMatchSeq ();
		
		if (start == startOriginal && stop == stopOriginal) return new int [] {oldStartQuery, oldEndQuery, oldStartMatch, oldEndMatch}; 
		
		int newStartQuery = oldStartQuery;
		int newEndQuery = oldEndQuery;
		int newStartMatch = oldStartMatch;
		int newEndMatch = oldEndMatch;
		
		int rebound = isReverse () ? (stopOriginal - stop) / 3 : (start - startOriginal) / 3; 
		int reduction = rebound - oldStartQuery;
		
		newStartQuery = Math.max (0, oldStartQuery - rebound);
		newEndQuery -= rebound;
		if (reduction >= 0) {
			newStartMatch = Math.max (0, oldStartMatch + reduction);
		} else {
			int d = Math.min (newStartQuery, oldStartMatch);
			int initMCompareSeq = oldStartMatch - d;
			int endMCompareSeq = oldStartMatch;
			if (endMCompareSeq > initMCompareSeq) {
				if (p.matchingSequence.length () - 1 < endMCompareSeq) {
					System.err.println ("WARNING: wrong endMCompareSeq " + initMCompareSeq + " " + endMCompareSeq);
					endMCompareSeq = p.matchingSequence.length () - 1;
				}
				if (endMCompareSeq < initMCompareSeq) {
					System.err.println ("WARNING: wrong initMCompareSeq " + initMCompareSeq + " " + endMCompareSeq);
					initMCompareSeq = endMCompareSeq;
				}
				String macthTail = OrfRecognizer.translateFromAminoToBase (p.matchingSequence.substring (initMCompareSeq, endMCompareSeq));
				String queryTail;
				if (isReverse ()) {
					int initQCompareSeq = stop + 3 - (newStartQuery - d) * 3;
					int endQCompareSeq = initQCompareSeq - d * 3;
					queryTail = OrfRecognizer.reverseComplementar (reco.wholeSequence.substring (endQCompareSeq, initQCompareSeq));
				}
				else {
					int initQCompareSeq = start + (newStartQuery - d) * 3;
					int endQCompareSeq = initQCompareSeq + d * 3;
					queryTail = reco.wholeSequence.substring (initQCompareSeq , endQCompareSeq);
				}
				if (getLevenshtein (queryTail, macthTail) > 80) {
	//				if (reduction < newStartMatch && reduction < newStartQuery) reduction = -Math.min (newStartMatch, newStartQuery); 
	//				else if (reduction < newStartMatch) reduction = -newStartMatch;
	//				else if (reduction < newStartQuery) reduction = -newStartQuery;
					newStartMatch -= d;
					newStartQuery -= d;
				}
			}
		}
		if (newEndQuery < newStartQuery) {
			newEndQuery = newStartQuery;
			newEndMatch = newStartMatch;
		}
//		if (newStartQuery < 0 || newEndQuery < 0 || newStartMatch < 0 || newEndMatch < 0) throw new RuntimeException (newStartQuery + " " + newEndQuery + " " + newStartMatch + " " + newEndMatch);
		return new int [] {newStartQuery, newEndQuery, newStartMatch, newEndMatch};
	}
	
	public float [] getMatchBox () {
		return getMatchBox (product);
	}
	
	public float [] getMatchBox (ProductProtein p) {
		if (!p.hasMatchingSequence ()) return new float [] {0, 0, 0, 0, 0, 0};
		
		int macthBound [] = getQueryMatchBounds (p);
		int startQuerySeq = macthBound [0];
		int endQuerySeq = macthBound [1];
		int startMatchSeq = macthBound [2];
//		int endMatchSeq = macthBound [3];
		
		int queryLength = length / 3;
		int matchLength = p.matchingSequence.length ();
		
		int boxBegin = Math.max (startMatchSeq, startQuerySeq);
		int boxEnd = boxBegin + (endQuerySeq - startQuerySeq);
		
		int startQ = boxBegin - startQuerySeq;
		int startM = boxBegin - startMatchSeq;
		
		int endQ = startQ + queryLength;
		int endM = startM + matchLength;
		
		float max = Math.max (endM, endQ);
		
		return new float [] {
			startQ /max, endQ / max,
			startM /max, endM / max,
			boxBegin / max, boxEnd / max
		};
	}
	
	public boolean isGene () {
		loadProduct ();
		return product.isGene ();
	}
	
	public boolean isFragment () {
		loadProduct ();
		return product.fragment;
	}
	
	public boolean inside (ORF orf) {
		return orf.start < start && orf.stop > stop;
	}
	
	public boolean sizesFit () {
		loadProduct ();
		if (!product.hasMatchingSequence ()) return false;
		float v = Math.abs (length / 3 - product.matchingSequence.length ()) ;
		return v < 2 || v / 	product.matchingSequence.length () < 0.001;
	}
	
	public float getSpecificity () {
		loadProduct ();
		if (!product.hasMatchingSequence ()) return 0;
		int macthBound [] = getQueryMatchBounds ();
		return Math.min (product.identity * (macthBound [3] - macthBound [2]) / product.matchingSequence.length (), 100);
	}
	
//	public float getAccuracy () {
//		if (!product.hasMatchingSequence ()) return 0;
//		int bb [] = getQueryMatchBounds ();
//		float a = (bb [1] - bb [0]);
//		float b = ((float) lenght / 3);
//		if (bb [1] < bb [0] || bb [3] < bb [2]) return 0;
//		float acc = (Math.min (a, b) / Math.max (a, b)) * (bb [3] - bb [2]) / (float) product.matchingSequence.length ();
//		return Math.min (getIdentity () * acc, 100);
//	}
//	
	public float getAccuracy () {
		if (product == null || !product.hasMatchingSequence ()) return 0;
		int bb [] = getQueryMatchBounds ();
		
		float lq = ((float) length / 3);
		int sq = bb [0];
		int eq = bb [1];
		float mq = eq - sq - sq - (eq - lq);
		if (mq < 0 ) return 0;
		
		float lm = (float) product.matchingSequence.length ();
		int sm = bb [2];
		int em = bb [3];
		float mm = em - sm - sm - (em - lm);
		
		if (mm < 0) return 0;
		
		if (bb [1] < bb [0] || bb [3] < bb [2]) return 0;
		float acc1 = (mq / lq) * (mm / lm) * getIdentity();
		float acc2 = ((eq - sq) / lq) * ((em -sm) / lm) * getIdentity();
//		float acc2 = 100;
		return Math.min (Math.min (acc1, 100), acc2);
	}
	
	private boolean loaded = false;

	public int overlaps;
	public void reloadProduct () {
		setLoaded (false);
		loadProduct ();
	}
	public void loadProduct () {
		if (isLoaded ()) return;
		setLoaded (true);
		File f = new File (CodonConfiguration.workspace + reco.dataDestDirS + "/orfs/" + id +  ".swiss");
		product = new ProductProtein (reco.dataDestDirS, id, f.exists () && f.length () != 0);
	}
	
	public boolean isUncharacterized () {
		loadProduct ();
		return product.isUncharacterized ();
	}
	public String getAminoSequence () {
		return OrfRecognizer.translateFromBaseToAmino (getSequence (), isReverse ());
	}
	public String getOriginalAmidoSequence () {
		return OrfRecognizer.translateFromBaseToAmino (getOriginalSequence (), isReverse ());
	}
	public String getSequence () {
		return reco.wholeSequence.substring (start, stop + 3);
	}
	public String getOriginalSequence () {
		return reco.wholeSequence.substring (startOriginal, stopOriginal + 3);
	}
	
	public String [] getTails () {
		loadProduct ();
		int startMatchSeq = getQueryMatchBounds (product) [2];
		String resp [] = new String [2];
		if (startMatchSeq > 2) {
			if (product.matchingSequence.length () < startMatchSeq) return resp;
			resp [1] = OrfRecognizer.translateFromAminoToBase (product.matchingSequence.substring (1, startMatchSeq));
			if (isReverse ()) {
				resp [0] = OrfRecognizer.reverseComplementar (reco.wholeSequence.substring (stop + 3, Math.min (stop + startMatchSeq * 3, reco.wholeSequence.length ())));
			}
			else {
				if (start - startMatchSeq * 3 + 3 < 0) return new String [] {null, null};
				resp [0] = reco.wholeSequence.substring (start - startMatchSeq * 3 + 3, start);
			}
		}
		return resp;
	}
	
	public void printTail (BufferedWriter writer) throws IOException {
		String [] tails = getTails ();
		if (tails [0] != null && tails [1].length () > 0) {
			LevenshteinResults res = levDDist.apply (tails [0], tails [1]);
			writer.write (id);
			writer.newLine ();
			writer.write (getIdentity () + "");
			writer.newLine ();
			writer.write (reco.wholeSequence.substring (Math.max (start - tails [0].length () - 50, 0), stop));
			writer.newLine ();
			writer.write (tails [0]);
			writer.newLine ();
			writer.write (tails [1]);
			writer.newLine ();
			
			int c = 0;
			String ss = tails [1];
			while (ss.contains (".")) {
				ss = ss.substring (ss.indexOf ('.') + 1);
				c ++;
			}
			
			int distance =  res.getDistance () - c;
			int substitutions = res.getSubstituteCount () - c;
			if (substitutions < 0) distance -= substitutions;
			writer.write ("D=" + distance);
			writer.write (" S=" + substitutions);
			writer.write (" I=" + res.getInsertCount ());
			writer.write (" D=" + res.getDeleteCount ());
			distance -= res.getDeleteCount ();  
			
			double lIndice = 100 - distance * 100. / (tails [0].length () - c);
			writer.write (" L=" + lIndice);
			writer.newLine ();
			String stExc = getIdentity () + ";" + lIndice;
			stExc = stExc.replace ('.', ',');
//			if (lIndice == 100) {
//				System.err.println ("(fragments) " + orf.id);
//			}
//			if (distance == 1 &&  queues [0].charAt (0) != queues [1].charAt (0) && queues [1].charAt (0) != '.') {
//				singleAlterationCount ++;
//				System.err.println (queues [0].charAt (0) + "*" + queues [1].substring (0, 3) + " " + orf.id);
//			}
//			if (distance == 3 && res.getInsertCount () == 3) {
//				System.err.println ("Codon insert " + orf.id);
//			}
//			System.err.println (stExc);
		}
	}
	
	public static int levCut = 500; 
	public static double getLevenshtein (String bases, String fromAminoToBase) {
		if (bases != null && fromAminoToBase.length () > 0) {
			if (bases.length () > levCut) {
				try {
					double lev1 = getLevenshtein (bases.substring (0, levCut), fromAminoToBase.substring (0, levCut)) * levCut;
					double lev2 = getLevenshtein (bases.substring (levCut), fromAminoToBase.substring (levCut)) * (bases.length () - levCut);
					return (lev1 + lev2) / bases.length ();
				}
				catch (StringIndexOutOfBoundsException e) {
					System.err.println ("Erro " + bases.length ());
					System.err.println ("fromAminoToBasek " + fromAminoToBase.length ());
				}
			}
			LevenshteinResults res = levDDist.apply (bases, fromAminoToBase);
			int c = 0;
			String ss = fromAminoToBase;
			while (ss.contains (".")) {
				ss = ss.substring (ss.indexOf ('.') + 1);
				c ++;
			}
			int distance = res.getDistance () - c;
			int substitutions = res.getSubstituteCount () - c;
			if (substitutions < 0) distance -= substitutions;
			distance -= res.getDeleteCount (); 
			return 100 - distance * 100. / (bases.length () - c);
		}
		return 0;
	}
	
	public double getTailLevenshtein () {
		String [] tails = getTails ();
		try {
			return getLevenshtein (tails [0], tails [1]);
		}
		catch (Exception e) {
			System.err.println (this);
			System.err.println (tails [0].length () + " " + tails [0]);
			System.err.println (tails [1].length () + " " + tails [1]);
			e.printStackTrace ();
			throw new RuntimeException ();
		}
	}
	
	public float getIdentity () {
		loadProduct ();
		return product.identity;
	}
	
	public void autoRebound () {
		List <Codon> stopStart = reco.startStopPorFrame [frame];
		loadProduct ();
		if (id.equals ("F1-001116-699-1815")) System.err.println ("Checking");
		if (product.identity < 80) return;
		int macthBound [] = getQueryMatchBounds (product);
		int startQuerySeq = macthBound [0];
//		int endQuerySeq = macthBound [1];
		int startMatchSeq = macthBound [2];
//		int endMatchSeq = macthBound [3];
		
		if (startMatchSeq < startQuerySeq) {
			if (!isReverse ()) {
				int ref = start + (startQuerySeq - startMatchSeq) * 3;
				int indexSS = reco.getCodonIndex (stopStart, ref + 2);
				if (indexSS == -1) return;
				Codon codon = stopStart.get (indexSS);
				if (codon.index > start && codon.index < ref + 3) {
//					int dd = (codon.index - start) / 3;
					start = codon.index;
//					product.endQuerySeq -= dd; 
//					product.startQuerySeq -= dd;
					length = stop - start;
				}
			}
			else {
				int ref = stop - (startQuerySeq - startMatchSeq) * 3;
				int indexSS = reco.getCodonIndex (stopStart, ref - 3);
				if (indexSS >= stopStart.size ()) return;
				Codon codon = stopStart.get (indexSS + 1);
				if (codon.index < stop && codon.index > ref - 3) {
//					int dd = (stop - codon.index) / 3 + 1;
					stop = codon.index;
					length = stop - start;
//					product.endQuerySeq -= dd;
//					product.startQuerySeq -= dd;
				}
			}
		}
		if (startMatchSeq + 1 > startQuerySeq) {
			if (getTailLevenshtein () > 80) {
				if (!isReverse ()) {
					int ref = start - (startMatchSeq - startQuerySeq) * 3 - 3;
					int index = reco.getCodonIndex (stopStart, start);
					if (index == -1) return;
					Codon codon = stopStart.get (index);
					while (codon.index >= ref && codon.type != OrfRecognizer.STOP) {
						index --;
						if (index == -1) break;
						codon = stopStart.get (index);
					}
					if (id.equals ("F1-001116-699-1815")) System.err.println (index);
					codon = stopStart.get (index + 1);
					if (codon.index < start) {
//						int dd = (start - codon.index) / 3;
//						product.setStartMatchSeq (startMatchSeq - dd);
//						product.setEndQuerySeq (startQuerySeq + dd);
						start = codon.index;
						length = stop - start;
					}
				}
				else {
					int refO = stop;
					int refN = stop + (startMatchSeq - startQuerySeq) * 3 + 3;
					int indexO = reco.getCodonIndex (stopStart, refO);
					int indexN = reco.getCodonIndex (stopStart, refN);
					if (indexN >= stopStart.size ()) return;
					Codon codon = stopStart.get (indexO);
					for (int i = indexO; i <= indexN; i++) {
						Codon test = stopStart.get (i);
						if (test.type == OrfRecognizer.STOP_R) break;
						codon = test;
					}
					if (codon.index > stop) {
//						int dd = (codon.index - stop) / 3;
//						product.setStartMatchSeq (startMatchSeq - dd);
//						product.setEndQuerySeq (startQuerySeq + dd);
						stop = codon.index;
						length = stop - start;
					}
				}
			}
		}
	}
	
	public void beheading (int size) {
		size = (size + 2) / 3;
		if (!isReverse()) stop -= size * 3;
		else start += size * 3;
		length = stop - start;
		if (length < CodonConfiguration.min_orf_length) {
			setRemoved (true);
			return;
		}
//		String s = getAmidoSequence ();
//		product.endMatchSeq -= size;
//		if (product.endQuerySeq > s.length ()) {
//			int diff = s.length () - product.endQuerySeq;  
//			product.endQuerySeq = s.length ();
//			product.endMatchSeq -= diff;
//		}
	}
	
	public int castrating (int size) {
		size = (size + 2) / 3;
		String s = getAminoSequence ();
		for (int i = size; i < s.length (); i ++) {
			if (s.charAt (i) == 'M' || s.charAt (i) == 'V') {
				if (!isReverse ()) start += i * 3;
				else stop -= i * 3;	
				length = stop - start;
//				product.startMatchSeq += i;
//				product.endQuerySeq -= i;
				if (length < CodonConfiguration.min_orf_length) setRemoved (true);
				return i;
			}
		}
		if (!isReverse ()) start = stop;
		else stop = start;	
		length = 0;
		setRemoved (true);
		return 0;
	}
	
	public boolean isBlasted () {
		loadProduct ();
		return product.blasted;
	}
	
	public boolean isFullBlasted () {
		File f = new File (CodonConfiguration.workspace + reco.dataDestDirS + "/full/" + id + ".unip");
		return f.exists ();
	}
	
	public boolean hasSwissUpdate () {
		File f = new File (CodonConfiguration.workspace + reco.dataDestDirS + "/orfs/" + id + ".swiss");
		return f.exists ();
	}
	
	public void saveAlterations () {
		String lg = "" + length;
		for (int k = lg.length (); k < 6; k ++) lg = '0' + lg; 
		String idNovo = "F" + frame + "-" + lg + "-" + start + "-" + stop;
		
		File s = new File (CodonConfiguration.workspace + reco.dataDestDirS + "/orfs/" + id + ".swiss");
		File d = new File (CodonConfiguration.workspace + reco.dataDestDirS + "/orfs/" + idNovo + ".swiss");
		if (!s.exists ()) {
			s = new File (CodonConfiguration.workspace + reco.dataDestDirS + "/orfs/" + id + ".unip");
			d = new File (CodonConfiguration.workspace + reco.dataDestDirS + "/orfs/" + idNovo + ".unip");
		}
		if (s.exists ()) {
			try {
				List <ProductProtein> products = ProductProtein.loadAll (s);
				for (int i = products.size () - 1; i > 0; i --) {
					ProductProtein p = products.get (i);
					try {
						int [] bounds = getQueryMatchBounds (p);
						p.setStartQuerySeq (bounds [0]);
						p.setEndQuerySeq (bounds [1]);
						p.setStartMatchSeq (bounds [2]);
						p.setEndMatchSeq (bounds [3]);
					}
					catch (Exception e) {
						products.remove (i);
					}
				}
				ProductProtein.saveAll (products, d);
				int [] bounds = getQueryMatchBounds ();
				product.setStartQuerySeq (bounds [0]);
				product.setEndQuerySeq (bounds [1]);
				product.setStartMatchSeq (bounds [2]);
				product.setEndMatchSeq (bounds [3]);
				startOriginal = start;
				stopOriginal = stop;
			} catch (Exception e) {
				System.err.println (id + " error while saving");
				e.printStackTrace ();
			}
		}
		id = idNovo;
	}

	public void reboundLeft (List <Codon> stopStart) {
		if (isReverse ()) {
			int index = reco.getCodonIndex (stopStart, stop) - 1;
			Codon codon = stopStart.get (index);
			if (codon.index > start) {
//				int dd = (codon.index - stop) / 3;
//				product.startMatchSeq -= dd;
//				product.endQuerySeq += dd;
//				if (product.startQuerySeq > 0) {
//					product.startQuerySeq -= product.startMatchSeq;
//					product.startMatchSeq = 0;
//				}
				stop = codon.index;
				length = stop - start;
			}
		}
		else {
			int index = reco.getCodonIndex (stopStart, start) - 1;
			if (index < 0) return;
			Codon codon = stopStart.get (index);
//			int dd = (start - codon.index) / 3;
//			product.startMatchSeq -= dd;
//			if (product.startMatchSeq < 0) {
//				product.startQuerySeq -= product.startMatchSeq;
//				product.startMatchSeq = 0;
//			}
//			product.endQuerySeq += dd;
			start = codon.index;
			length = stop - start;
		}
	}

	public void reboundRight (List <Codon> stopStart) {
		if (!isReverse ()) {
			int index = reco.getCodonIndex (stopStart, start) + 1;
			Codon codon = stopStart.get (index);
//			int dd = (start - codon.index) / 3;
			if (codon.index > stop) return;
//			product.startMatchSeq -= dd;
//			if (product.startQuerySeq > 0) {
//				product.startQuerySeq -= product.startMatchSeq;
//				product.startMatchSeq = 0;
//			}
//			product.endQuerySeq += dd;
			start = codon.index;
			length = stop - start;
		}
		else {
			int index = reco.getCodonIndex (stopStart, stop) + 1;
			Codon codon = stopStart.get (index);
//				int dd = (codon.index - stop) / 3;
//				product.startMatchSeq -= dd;
//				product.endQuerySeq += dd;
//				if (product.startMatchSeq < 0) {
//					product.startQuerySeq -= product.startMatchSeq;
//					product.startMatchSeq = 0;
//				}
				stop = codon.index;
				length = stop - start;
		}
	}

	public boolean isRemoved () {
		return removed;
	}

	public void setRemoved (boolean removed) {
		if (isExcluded ()) removed = true;
		if (isForced ()) removed = false;
		this.removed = removed;
	}

	public boolean isExcluded () {
		return excluded;
	}

	public void setExcluded (boolean excluded) {
		this.excluded = excluded;
		if (excluded) {
			forced = false;
			removed = true;
		}
		else {
			removed = false;
		}
	}
	
	public boolean isForced () {
		return forced;
	}

	public void setForced (boolean forced) {
		this.forced = forced;
		if (forced) {
			this.removed = false;
			this.excluded = false;
		}
	}
	
	public String toString () {
		return this.id + "�" + this.startOriginal + "�" + this.stopOriginal + "�" + this.start + "�" + this.stop + "�" + this.removed + "�" + this.forced + "�" + this.excluded + "�" + (this.product == null ? "null" : this.product.toString ());
	}
	
	public static ORF createFromString (OrfRecognizer reco, int frame, String line, int version) {
		String args [] = line.split ("�");
		if (version == 1) {
			String tmps [] = new String [args.length + 2];
			for (int i = 0; i < 3; i++) tmps [i] = args [i];
			for (int i = 1; i < args.length; i++) tmps [i + 2] = args [i];
			args = tmps;
		}
		ORF orf = new ORF (reco, Integer.parseInt (args [1]), Integer.parseInt (args [2]), frame);
		orf.start = Integer.parseInt (args [3]);
		orf.stop = Integer.parseInt (args [4]);
		orf.id = args [0];
		orf.setForced (Boolean.parseBoolean (args [6]));
		orf.setExcluded (Boolean.parseBoolean (args [7]));
		orf.setRemoved (Boolean.parseBoolean (args [5]));
		if (args.length > 9) {
			orf.product = ProductProtein.createFromStrings (args);
			orf.setLoaded (true);
		}
		return orf;
	}

	public boolean isLoaded () {
		return loaded;
	}

	public void setLoaded (boolean loaded) {
		this.loaded = loaded;
	}

	public boolean isReverse () {
		return frame > 2;
	}
	
	public boolean isForward () {
		return frame <= 2;
	}

	public void selectBetterEntryThatFit () {
		File f = new File (CodonConfiguration.workspace + reco.dataDestDirS + "/orfs/" + id + (hasSwissUpdate () ? ".swiss" : ".unip"));
		if (!f.exists ()) return;
		try {
			List <ProductProtein> products = ProductProtein.loadAll (f);
			if (products.size () == 0) return;
			float max = 0;
			for (ProductProtein p: products) {
				ORF orf1 = new ORF (reco, startOriginal, stopOriginal, frame);
				orf1.product = p;
				orf1.setLoaded (true);
				float acc = orf1.getAccuracy ();
				if (max < acc) max = acc;
			}
			final float maxAcc = max;
			Collections.sort (products, new Comparator <ProductProtein> () {
				public int compare (ProductProtein p1, ProductProtein p2) {
					ORF orf1 = new ORF (reco, startOriginal, stopOriginal, frame);
					ORF orf2 = new ORF (reco, startOriginal, stopOriginal, frame);
					orf1.product = p1;
					orf1.setLoaded (true);
					orf2.product = p2;
					orf2.setLoaded (true);
					
					float acc1 = orf1.getAccuracy ();
					float acc2 = orf2.getAccuracy ();
					
//					if (acc1 < 20 && acc2 < 20) return acc1 > acc2 ? -1 : (acc1 == acc2 ? 0 : 1);
					if (acc1 < maxAcc - CodonConfiguration.tolerance && acc2 < maxAcc - CodonConfiguration.tolerance) return acc1 > acc2 ? -1 : (acc1 == acc2 ? 0 : 1);
					if (acc1 >= maxAcc - CodonConfiguration.tolerance && acc2 < maxAcc - CodonConfiguration.tolerance) return -1;
					if (acc1 < maxAcc - CodonConfiguration.tolerance && acc2 >= maxAcc - CodonConfiguration.tolerance) return 1;
					
					if (!p1.fragment && p2.fragment) return -1;
					if (p1.fragment && !p2.fragment) return 1;
					
					if (p1.isGene () && !p2.isGene ()) return -1;
					if (!p1.isGene () && p2.isGene ()) return 1;
					
					if (!p1.isUncharacterized () && p2.isUncharacterized ()) return -1;
					if (p1.isUncharacterized () && !p2.isUncharacterized ()) return 1;
					if (orf1.length > orf2.length) return -1;
					
					if (orf1.length < orf2.length) return 1;
					if (acc1 == acc2) return 0;
					if (acc1 > acc2) return -1;
					if (acc1 < acc2) return 1;
					return 0;
				}
			});
//			ProductProtein.saveAll (products, f);
			product = products.get (0);
		} catch (IOException e) {
			e.printStackTrace ();
		}
	}
	

	public void selectBetterEntryFulfill (int freespan) {
		File f = new File (CodonConfiguration.workspace + reco.dataDestDirS + "/orfs/" + id + (hasSwissUpdate () ? ".swiss" : ".unip"));
		if (!f.exists ()) return;
		try {
			List <ProductProtein> products = ProductProtein.loadAll (f);
			products.sort (new Comparator <ProductProtein> () {
				public int compare (ProductProtein o1, ProductProtein o2) {
					if (o1.matchingSequence.length () > o2.matchingSequence.length ()) return -1;
					if (o1.matchingSequence.length () < o2.matchingSequence.length ()) return 1;
					return 0;
				}
			});
			if (products.size () == 0) return;
			float max = 0;
			for (ProductProtein p: products) {
				ORF orf1 = new ORF (reco, startOriginal, stopOriginal, frame);
				orf1.product = p;
				orf1.start = start;
				orf1.stop = stop;
				orf1.setLoaded (true);
				orf1.autoRebound ();
				float acc = orf1.getAccuracy ();
				if (max < acc) max = acc;
			}
			final float maxAcc = max;
			
			float selAccuracy = getAccuracy ();
			ORF selOrf = this;
			product = products.get (0);
			autoRebound ();
			for (ProductProtein p: products) {
				if (product == p) continue;
				
				ORF orf1 = new ORF (reco, startOriginal, stopOriginal, frame);
				orf1.start = start;
				orf1.stop = stop;
				orf1.product = p;
				orf1.setLoaded (true);
				orf1.autoRebound ();
				
				if (	product.isGene () && 
						selAccuracy > maxAcc - CodonConfiguration.tolerance && 
						!product.fragment && 
						selOrf.length < freespan && orf1.length < length) {
					break;
				}
				
				float acc1 = orf1.getAccuracy ();
				if (acc1 < maxAcc - CodonConfiguration.tolerance) {
					continue;
				}
				if (acc1 >= maxAcc - CodonConfiguration.tolerance && selAccuracy < maxAcc - CodonConfiguration.tolerance) {
					product = p;
					selAccuracy = orf1.getAccuracy ();
					selOrf = new ORF (reco, startOriginal, stopOriginal, frame);
					selOrf.start = start;
					selOrf.stop = stop;
					selOrf.product = p;
					selOrf.setLoaded (true);
					selOrf.autoRebound ();
					continue;
				}
				if (!p.fragment && product.fragment) {
					product = p;
					selAccuracy = orf1.getAccuracy ();
					selOrf = new ORF (reco, startOriginal, stopOriginal, frame);
					selOrf.start = start;
					selOrf.stop = stop;
					selOrf.product = p;
					selOrf.setLoaded (true);
					selOrf.autoRebound ();
					continue;
				}
				if (p.isGene () && !product.isGene ()) {
					product = p;
					selAccuracy = orf1.getAccuracy ();
					selOrf = new ORF (reco, startOriginal, stopOriginal, frame);
					selOrf.start = start;
					selOrf.stop = stop;
					selOrf.product = p;
					selOrf.setLoaded (true);
					selOrf.autoRebound ();
					continue;
				}
				if (!p.isUncharacterized () && product.isUncharacterized ()) {
					product = p;
					selAccuracy = orf1.getAccuracy ();
					selOrf = new ORF (reco, startOriginal, stopOriginal, frame);
					selOrf.start = start;
					selOrf.stop = stop;
					selOrf.product = p;
					selOrf.setLoaded (true);
					selOrf.autoRebound ();
					continue;
				}
				if (orf1.length < freespan && selOrf.length < freespan) {
					if (orf1.length > selOrf.length) {
						product = p;
						selAccuracy = orf1.getAccuracy ();
						selOrf = new ORF (reco, startOriginal, stopOriginal, frame);
						selOrf.start = start;
						selOrf.stop = stop;
						selOrf.product = p;
						selOrf.setLoaded (true);
						selOrf.autoRebound ();
						continue;
					}
				}
				if (orf1.length < freespan && selOrf.length > freespan) {
					product = p;
					selAccuracy = orf1.getAccuracy ();
					selOrf = new ORF (reco, startOriginal, stopOriginal, frame);
					selOrf.start = start;
					selOrf.stop = stop;
					selOrf.product = p;
					selOrf.setLoaded (true);
					selOrf.autoRebound ();
					continue;
				}
				if (orf1.length > freespan && selOrf.length > freespan) {
					if (orf1.length < selOrf.length) {
						product = p;
						selAccuracy = orf1.getAccuracy ();
						selOrf = new ORF (reco, startOriginal, stopOriginal, frame);
						selOrf.start = start;
						selOrf.stop = stop;
						selOrf.product = p;
						selOrf.setLoaded (true);
						selOrf.autoRebound ();
						continue;
					}
				}
				if (acc1 > selAccuracy) {
					product = p;
					selAccuracy = orf1.getAccuracy ();
					selOrf = new ORF (reco, startOriginal, stopOriginal, frame);
					selOrf.start = start;
					selOrf.stop = stop;
					selOrf.product = p;
					selOrf.setLoaded (true);
					selOrf.autoRebound ();
					continue;
				}
			}
			
//			Collections.sort (products, new Comparator <ProductProtein> () {
//				public int compare (ProductProtein p1, ProductProtein p2) {
//					ORF orf1 = new ORF (reco, startOriginal, stopOriginal, frame);
//					orf1.start = start;
//					orf1.stop = stop;
//					orf1.product = p1;
//					orf1.setLoaded (true);
//					orf1.autoRebound ();
//					
//					ORF orf2 = new ORF (reco, startOriginal, stopOriginal, frame);
//					orf2.start = start;
//					orf2.stop = stop;
//					orf2.product = p2;
//					orf2.setLoaded (true);
//					orf2.autoRebound ();
//					
//					float acc1 = orf1.getAccuracy ();
//					float acc2 = orf2.getAccuracy ();
//					
//					if (acc1 < maxAcc - CodonConfiguration.tolerance && acc2 < maxAcc - CodonConfiguration.tolerance) return acc1 > acc2 ? -1 : (acc1 == acc2 ? 0 : 1);
//					if (acc1 >= maxAcc - CodonConfiguration.tolerance && acc2 < maxAcc - CodonConfiguration.tolerance) return -1;
//					if (acc1 < maxAcc - CodonConfiguration.tolerance && acc2 >= maxAcc - CodonConfiguration.tolerance) return 1;
//					
//					if (!p1.fragment && p2.fragment) return -1;
//					if (p1.fragment && !p2.fragment) return 1;
//					
//					if (p1.isGene () && !p2.isGene ()) return -1;
//					if (!p1.isGene () && p2.isGene ()) return 1;
//					
//					if (!p1.isUncharacterized () && p2.isUncharacterized ()) return -1;
//					if (p1.isUncharacterized () && !p2.isUncharacterized ()) return 1;
//					
//					if (orf1.length < freespan && orf2.length < freespan) {
//						if (orf1.length > orf2.length) return -1;
//						if (orf1.length < orf2.length) return 1;
//					}
//					if (orf1.length < freespan && orf2.length > freespan) return -1;
//					if (orf1.length > freespan && orf2.length < freespan) return 1;
//					if (orf1.length > freespan && orf2.length > freespan) {
//						if (orf1.length > orf2.length) return 1;
//						if (orf1.length < orf2.length) return -1;
//					}
//					
//					if (acc1 == acc2) return 0;
//					if (acc1 > acc2) return -1;
//					if (acc1 < acc2) return 1;
//					return 0;
//				}
//			});
//			product = products.get (0);
			autoRebound ();
		} catch (IOException e) {
			e.printStackTrace ();
		}
	}
	
	public void selectBetterEntryMaximizingIntergenicRegion () {
		File f = new File (CodonConfiguration.workspace + reco.dataDestDirS + "/orfs/" + id + (hasSwissUpdate () ? ".swiss" : ".unip"));
		if (!f.exists ()) return;
		try {
			List <ProductProtein> products = ProductProtein.loadAll (f);
			if (products.size () == 0) return;
			float max = 0;
			for (ProductProtein p: products) {
				ORF orf1 = new ORF (reco, startOriginal, stopOriginal, frame);
				orf1.product = p;
				orf1.start = start;
				orf1.stop = stop;
				orf1.setLoaded (true);
				orf1.autoRebound ();
				float acc = orf1.getAccuracy ();
				if (max < acc) max = acc;
			}
			final float maxAcc = max;
			Collections.sort (products, new Comparator <ProductProtein> () {
				public int compare (ProductProtein p1, ProductProtein p2) {
					ORF orf1 = new ORF (reco, startOriginal, stopOriginal, frame);
					orf1.start = start;
					orf1.stop = stop;
					orf1.product = p1;
					orf1.setLoaded (true);
					orf1.autoRebound ();
					
					ORF orf2 = new ORF (reco, startOriginal, stopOriginal, frame);
					orf2.start = start;
					orf2.stop = stop;
					orf2.product = p2;
					orf2.setLoaded (true);
					orf2.autoRebound ();
					
					float acc1 = orf1.getAccuracy ();
					float acc2 = orf2.getAccuracy ();
					
					if (acc1 < maxAcc - CodonConfiguration.tolerance && acc2 < maxAcc - CodonConfiguration.tolerance) return acc1 > acc2 ? -1 : (acc1 == acc2 ? 0 : 1);
					if (acc1 >= maxAcc - CodonConfiguration.tolerance && acc2 < maxAcc - CodonConfiguration.tolerance) return -1;
					if (acc1 < maxAcc - CodonConfiguration.tolerance && acc2 >= maxAcc - CodonConfiguration.tolerance) return 1;
					
					if (!p1.fragment && p2.fragment) return -1;
					if (p1.fragment && !p2.fragment) return 1;
					
					if (p1.isGene () && !p2.isGene ()) return -1;
					if (!p1.isGene () && p2.isGene ()) return 1;
					
					if (!p1.isUncharacterized () && p2.isUncharacterized ()) return -1;
					if (p1.isUncharacterized () && !p2.isUncharacterized ()) return 1;
					if (orf1.length < orf2.length) return -1;
					
					if (orf1.length > orf2.length) return 1;
					if (acc1 == acc2) return 0;
					if (acc1 > acc2) return -1;
					if (acc1 < acc2) return 1;
					return 0;
				}
			});
			product = products.get (0);
			autoRebound ();
		} catch (IOException e) {
			e.printStackTrace ();
		}
	}
	
//	public ORF () {}
//	public static void main(String[] args) {
//		String a = "GTTCAGATTCAGGAAGAACTTCGCGATCTGCCTGGCCCGCGTCTGGCGATGACTCAACCCATCGAAATGCGGATGAATGAAATGATCTCCGGGGTCCGCTCAGATGTAGCTGCGATTCTCTACGGAGACGATCTCGATCTGATGGTCTCCAAAGCCTCTGAAATTGAGCAGGTGCTCAAATCAATTGAAGGGGCTGCCGACGTGAAAGTGGAGCAGGTCACCGGTCAGCCGGTGCTGGAAATCAACATCGACCAGGATGAAATCGCCCGTTACGGCATTCCGGCACGGACCGTGCTCGACCTGGTGGAATCATTGGGCAGCAAACATGTCGGAGAGGTCTATGAAGGACAGCTGCGATTTCCCCTGATTATCCGGCTACCCGAGCAGGCACGCACGGATCCCGAGGCGATCGGCTCGATCCTGGTAGCCACTCCTCAGGGAGAGCAGATTCCCCTGTCGCGACTCGCTGACATTGAGATCGTGGAAGGCCCGAATACCATCAAGCGGGAATGGTATCAGAGGCGGATTACCATCGAAGCAAATGTCCGGGGTCGCGATATGGGGAGCTTTGTCGCTGAAGCGCAGCGTAAAGTCGACAAGCAGATAACCCTGCCAGCCGGACGTTACCACATCGAATGGGGTGGACAGTTCGAGAACCTGCAGCGTGCCCAGGCCCGCCTGATGATTGTGGTGCCGGTCGCGATGCTACTGATTTTCGGCTTGCTGTATATGACTTATAATAACCTGATCGATTCAATTCGCGTCTTTACTGGCGTCCCTTTCGCCTGGGTCGGGGGAATTTTCGCACTCTGGATACGGGAAATGCCGTTTTCGATTTCAGCGGCCGTGGGTTTTATCGCCCTCTCCGGTGTCGCGGTACTGGATGATATGCTTCTGGTCTCAACCATCCGTCGATTGCGCCGCAGAGGCAGTGCGCTAAACGAAGCTGTCGAAGAGGCAGCTATGACCCGTTTGCGTCCGATTCTGATGACTACGCTGGTGGCCAGCCTCGGCTTCTTTCCAATGGCCTTCAATACCGGTATGGGAGCAGAAGTGCAACGCCCCCTGGCCACAGTCGTGATTGGTGGCGTCTGCAGTGCCACGATCATGAGCCTGCTCGTTTTACGTGTGCTGTACGTGGTCTTTAATTTACCCGCAAGTCGCAACGACGATGAGGAAGACGATGACGATACCACTAAACCAACACAACCGAACAATCCGGATTCCAGGCATGTGCCTGAAACGGCTTCGGTTTAATATGACTGGTCAATGTACCAAAAAAGGAGAAAAAATGACTAAACAATCATGGACCCGGTGTCTGACAATCCCTGCGATTGTGTGTGGACTTCTGGTGATCACGGGCTGTGCGGAGAGTGGAAATGACGTAGCCCCCGTGGCCAAAGCTGAAGAAGACCATCATGAGCATGGCGAAGAAGGGCACGATCATGAACACGAACACGGAGAAGAGGGGCATAACCTGCACGGTTTCTGGTGTGCTGAGCATGGCGTTCCGGAGGAAATCTGCGCCCAGTGTAGTACCAAACTCGCCGCTGATTTCCAGAAGAAGGGAGACTGGTGTGAAGAGCATAAACGTCCTGATTCGCAGTGCTTCATCCATCATCCCGAACTGGAAGAAAAATTCATTGCCCAGTACGAGGCCAAGTTCGGTGAAAAACCGCCCGCACGTAAGGATAAATAAGATCCGCGAGCAGACAAAAGCGGAAACACAGTCGATTAGCAGCCAGGTAACCCGCCCGAATTCAAGCGGGTTGCCTGTTCTGTTTTACGCTGTCGCACTGCTGATCTGCGCATTGATTCCTGGCTGCCAGAGTGCCGCTTCAGATTCATCTCATTCCAGTTTGCCTGAGGTTACAGATTCGAATTTTCAAAAGTCGGTCCTGGAGGCGGACCAGCCAGTGCTGGTGGAATTCTGGGCTCCCTGGTGCCGCCCCTGTATCGAGATGATCCCCCTCCTGGAAGAGGCGTCAGAACAATTTGCGGGACGCGTGAAGATTCTCAGAATGCGCATCGACGAGAACCCGGCTACTGCTGCCAAATATGAAATCGACGCGCCCCCCGCGTTCCTTCTCTTTAACGAGGGGAAGGTCTTCAAGCGCCGACTCGGGAAACAAACACAAGCACAACTGACGAGCCTCATAGAGAATGTCTGCGCAGAAACTGCCTCGACTCGTCGGGTACAAATAGACTGTCGCCCGTTAAGTCTCGTAGCCAATTGATCATTGAAAACAAAATCTGGTGAAAGAGCTGACGTGGGTCGGGTCACATCTGCCGTGAATGACTATTGGTTGTGATATTTTAGGAGCGAACGCGCGTACCCTGGGCGTGAAGAAATGACTTGCTTTTCGATGTGATCTACGCCGCAATCCGCTCGTACGACTTGATCAGTCCTCCGAGATGTTCCTCGCAGTGGATCTCATCGAGCCTGATTGTCTCGAATTCCGGCGGTGGCTCAGCACAACAGGGCGGCAGATATTCTCGTGAGCTGTGCGCTCGATGTTTGTTGTAATAATCAACGAATTCCGAAACAAGATAATCGAGATGCACTTTGCCGAAAGGCTGCGTACGTTCGTGAAGCTTTGTATGAATTTCACCTGGTATACGATCCCGCAGAATCCGGTTCTCATTTTTCAAATAGAGGACGTATTTGGCCAACAGGCTGTCGGAAGCAGTCGCAATCAGTGTCAGCAGAGGATGAAAAAGTTGGATCATCGTGGTAAAATCAGTCTCGGTGGTAACATGCTGTAAAATAAGGAGATGTCACTTAATTTGTACCCCGCTGAGACACAGTTCACATTCGATGGACAGTCTCTTTCAGGGCTTTGCGGAACCTTCACCGACTGCAGGGCCACAGAGTTGACGCAGTTCACCTTATATAATGGTGTCGCGGCTGCTATTAAGCTACTCCAGAAAAATTCGAAAGCGGTGAGTTACCTCTAAATCACCGTAGAGGGCTCGCCCTTTTGTAAGCTCCTTTGATCTCATGACTTCCGGAATCAGCCCGCGCAATAGTCGTTACGAAGAGGGTTTTGTACGGGACGGGGAGCTCTACTTCTAATTGTTGAATGGGAAACAATCGAATTGCCTCCGGTTCGTGATTTTCGATGGCTTGGCCATACTCGAAACCATCGACTTTGGCGATATATTTACCTCCTGTGACACCTGCTTCCTCAAGTGTCGGTATAGGCCATTTTCAAAATCAAGGCTCAGTCGGGCCCCCTCATTACCTGCTACTCCATCCTCGAATTTCCTGTTCGATTTTGAACCAGATGCCGAGATTACTTTTCAGTCACTCCCGTAGCAGTAGTTTTCTCGTCCGGTCTGGCGGGAACAACGCGAATCACGCGACCGCTCTGGTGCTTGTCGCCATCCATTTCCACCGTGCGAACTCGCAGCAGATGCGTGCCCGGCTTCAAACCGGCGGGCAGTTTCGCATGCCAGAGATGTGTGGACTTTTTGGCATTCGGCATACTACGGTACTTCTTATCTTTCACACTGTTCTCAGTCTCAGACAGTTTCTTGAAGCCCGGATCAATGCCTGGCCGAAATTCCATCGGAACCCAGGTTCCAACCGCTCCGATCAGCATTGCAACTTTGGAACGTTCCGAACCATTGAACACGTTGACAATCACATCTGTCTCGACCAACTGGTTGGCAGCGACTTCTTCCGGAACCATGATATTCATCTGATACTTCTCAGAGCGGCCTGCCGCCCGAAAGTCCAGATTGTATTCCTTGCCATCAAAACGAATGATGGAATAACCATTCGGCGCCCCGTCGGCCATCACGGTATGCGGAATGCCCCGTTCGTCAGGAGTGCCAGACCACCAGCTGCCACAAACCGTTACGTTGATGATGTGATGATGCGGTTTCGGGCCACGCCAGCCATCTGCTTTGGTGATGAAACGGTGCTCGTGATGATGCGTGTGACCGGAAATCGACATGCAAAACGGACGTTTTTCAATCAACCGGTACAACTCCTGTCGATCGTTGACGTTCACCAGCGGGATATGCATCATCAGCACAACCAGCTGATCTTCGGGAATCTGTTTCAGGTCGTTCCTGATAAACTCCATCTGCTCTTTACCCAGGCCGCCTTCGTATTTTCCTTTCTTGTCTTTTTCTGAAACAGTCCATTCGACATCATCGAGGGTAATAAAATGCACGGTCCCATAATCGAATGAGTAATACGCCGGCCCATACTCACGTTCGAAAGTTTCATCACTCAGTTTGTCATTCGGGGCATCAAAGTTGATGTCGTGATTGCCGATCACGTTGTACCAGGGAATGCCCAGCAGCGCGATGCCGCGCGCCTGCGATTGGAACAGTGAAAGATCGTTAAATAGTATATCTCCCAGAGTCACGCCGAACGAAGCGTCGGTTCCAATCAGATCTTCGATGACATCGTGGGCGATATAATCAATTTCCTTTTGGTCACGTGGCTGCGGGTCGCCAAAAAAGATCGCCCGGAATTCAGCAGGCTCTTTCTGGGGATACAGAGGAAAATCCACCGATGCAGGCAGCGGACCGGTCGGCTTGACCCCCGCGTACTTTGATTTTGGTGAGCCGTTCGGCTTGTGAATATAATAGAATTCCGGCGTCAGATTCTTGCTCAGCGGCGTCCGCCATCCTTGTGGTTTGATTACAAACAGAATCGTATCATCCGTAACCGGCAGTTGATACAAACCATTTTCGTCCGTGCAGACAACCTCGCGACCATTGGAAACGCGAACACAAGCCAGGCCCTGTTCCCCGGCGTCCCGCTGGCGATTGTTGTTCGCATCATAAAACACATAACCGGTGGCTGTCTGCGACTCGTTTGCTGTTTTTGATTCCGCCCGTCCGTAGGTTTGACCACTCGATACCAGCAGTAAGACCACAAGCAATTTTGTCAGTGATGATCGACTCATTATTTTCCCTTATCAGTAAAAATTTGCTCTAATTATGTTATGTTCCGTAGTATCACGGTAAAACCTATAGACCATTATTTATGTGAAGTTTCCGGCAATACTATTCAAACCGTAAGTAGAAATCCTTACTAATTTCCTCATATATTGAGTATTATTACAATAACATCATACAAAAAGCAGAAAGGGGTTTCAGGGAAATACAGGAGCCTATCAAAAGCGACAAGCCGCGTCCGCAATACGAGCAGAAAAATCTCGAAAACGACGGATGGAATAACGTGGCCACTCGCCACGCTAAAAGTGAACCAAATGCAAGACAGCGCAAGAAATTCGGGTGCGAGGATCACCAGCTGTAAATTCTGTCGGCTTCAGTCGTCCGCACGATCTGACAGGATTTTCCCCGTCGACGCCTCAATAGTCAGCTCCCGCTCAGAACCGGCATGCAGCGTCTCTACTTCCCATTGACCATGCTCATAATCCATTTCCAGGATAGGCTGATAACCCGCCTGCTCTATCTGGCGGGCGATCTGGGAAGCGGCCAGGTCTTTGATACTCAGTTTCTGATGATGGAAATCCGGATGCTCGCTAATAATCGCCAGCGTCTGCGGGTTGACGTGGACCTCAATCGGCTTGTTCTTGCGAAAGGTTTCGATTTCCCAGACACCATCGTCAAACGCAGCTTCCACAATAGGTGTCAGCCCGAATGGTACTGTCATTTTCCGATTTGTGCGTACGTTGGTCTATGACAGAATCGAATCTTTTGTAAAATCACTTCTGTCTGGAGTTGGTCTGGCTTTTTAACCGCTCCCGGCGATAATGTTTTCAACCCTCCTTAGGTCTTCCGGCAAGAAATCACCAGGGAGGGTTTTTTCATGGCCGGTAATCTATCTGAATCGTCGTCAGAAGGCAGTAGCAGTTTTCAACTGATTCTCGATTCCTTCCTGTCATCGGCGGGACTACCTTTTTCTAACCTGCTCTCTGCCGAGCGCATTGCGCGGATCTTTGCACGGCATAACGGGTTGTTCGGAACTCACGGCGTTTATTCTACCGCCGTCATGGTCTGGTCTTTTCTGGGACAGGTTCTACGGGACGGAAAAGAAGCTTCCTGCCAATCTGCCGTGGCACGTGTCATGACCTACTGCCAGTTGACCGGCAAGACAACGCCGACCGCTGATACGCGTAACTATTGCCGGGCACGGGCAAAACTCTCTGTGCCGGCATTGCGTGAACTCAGCGGTGAAATTGCCAACGAACTCGAAGAAGCAGCCGAGGAGCACTGGCTCTGGAAAGGCAGGCATGCCAAACTGATCGACGGCTTCACTTTCACCATGCCGGACACAGAACTGAATCAGGATGAGTTCCCACAACAGTCTGCCCAATCTCCGGGCTGTGGTCTGCCGATCGCCCGCGCCGTGGCGGTGCTGTCACTGGCAACAGCCTGCGTGCTGGATGTCGCCATCGGTCCCTATAAGGGAAAGCAGACTGGTGAGACAGCTCTGCTGCGCAAACTGTTCGCAGCGCTGAAGCCCGGTGATGTTGCCATCTTCGATCGTTACTACTGCTCATTTATGATGATCGCTTCGCTACTCAACCAGGGAACCGATGTCTGTGCGCGACTGCATCACAAACGTCGTGCTGACTTTCGACGTGGCAGACGACTGGGCAAATACGATCATCTCGTTTCCTGGACCAAACCCATACAATGTCCCACCTGGATGGATCGGGTGACTTTCGAGCAGATCCCGGAAACAATACTGCTGCGGGAAGTACGTTACAATCTGGTGGAGCGAGGACGGAGGACTCAATCGATGACCGTCGTTACCACCCTCATCGACGAGGCGTTGACCAAAGAGGAAATTGCAGAACTATATGGGTTTCGCTGGAATGCAGAACTCGATCTGCGCTCGATCAAAGACACACTCAACCTGGGACACCTGCGCTGCAAGTCGCCGGCAATGATTCGCTGTGAACTCTGGACCACACTGCTGGGCTACAACCTGATTCGCACCACCGCTGCCGGTGCGGCGGTGTTACATGAGAAACGTCCCCGACAGATCAGCTTCACGGCAACGTGCCAGTTGATTCTATCGGCCTGGATGAACACGGCCTGCGGCAAATTGGAGCAGCAGAAATTAAAAGTACTCTGTCAAACCCTGGCAGATCAGATTGCCGTTTGCGAAGTGGGCAACCGGCCTGGCCGAATCGAACCCCGCGTTATCAAACGCAGACGCGGCACCTATCCACTGATGCAAAAACCCAGACAGGCTCTCAGAGACCATTTACGCAACAACACCACTTGAAACATCTTAGGACTCGACAGTACCATTCGTGTCAGCCCCTGTTGTTCCAGCTTTGATAAAACCTGTGAAAGCGGACTGGCCTGGCTCGGCGGGACATCAGCTTCTGCTACGGTTACTACTACGCTCACAAGACTTCCCAATAATAATACGGTACGCATGATTGACTCCTTAGTCAGTAAATATTTTGGAAAAGGATTTCCGTTCTTACAGTCAGATACTATCAAGCGTGCCTGAAATCAGGATGAACCAAAGATAAAATCCACTTCATGTTTGAAAAATAAATCAGCTTGCATTACCGGGCTTTATCAAAGCTGCCACCTAATATGAAATAGACCGGGAAACGGAAAAAATATTTGACAGCCATAATATGTCACATATCGTATGTCTGGATGTGAAGTCAACTATGACACCAGTCACACGACTTTGTGTCCGATTGATTCTTATCTTCCCGGTGACATCAGCCACCGGAGTCCTTGATTCTATCTGATGTAAGGTGAAGTCGAATTGTTCTTCAGTGTATCGTGTTCGTTTCATTTCCCGCTCCTTTCTCAAAATATAAGCTAACAAAACCTCTCACATTTCACATGGATCAAGAAACGGGGTGGGAAGATCACTCTATTTTCTGACACAGAATTAATCAGAACCAAAGCTCATTCCCGGCCTGGAACAGGTTGACACGCGGGCTTACCATGATGACTATCCGTCAGAACGGTGAAACACGACTGATTACATCAGCCGGTCATCGTTGCGTTTGAGTTCGTTAAAGTCTTCGGTTTCCGGTAAAAGGACGGCGAGCTTCTCACTTCTGTTAAGGCGAATCACGCTACTGGTCAGTATATCTTTCAGATAGTAACTGACCTCGACCGCTTCCCAATCCATGCGAAAGATGGAGTCGGTTATTCCCAGTGAATTACAGACATCGATTGTGGATTTTTGGCGAGAGATCAGTGTCAGTCGGCGGCCCAGGTTCTCGCGGCGGAGTAATTGATCGAACTGGAGTGGCAAGGATTCGTCAAGCAGCAGAACAACATCGGAAATCCCCTGCCCGTTTGGCTTATCGAAACGAGACAACACTACTCCCAGTTCGGCCGTAGAGGCATAGCGATGAATTTCCACATGTGAGGCACCCCAAGATGACCATTGCTTGATCGTATCATTCAGGTCCTGCATATTGGCGGATGAGCCGCTATGCTGTTTTTTAAATCGGGTTCGAGTGTCTCTTATATCGGATCGCTTTGTGATCGAGATCAGATCATTGATAAAACTGGTAAGAGTTACCATGATAAATTTCCTTTTTCTGGATGGCTATTCATATTATTGATTAAATGAGTGATTAATCTGGTAAGTGGAGTCTGGAACTCGTCACGATGAAGAGGATGCTGTGATTCGGTCACAAGACTTAGTTTGCTGCTGTTCCACGAAACCCTTTTTATTGTTTGTTTTCTGTTCAATTCGCTGCACCAGATACGGAATGGCGTCTTGAATTGTCAGACAGGTGTGATACTCGGGCCTGATTATTCTGAAGCGTTGAAGGTAAGGGATCATCTAGATTTCCTTTGATCAGTAGACTATGAAATTTGCCGGTGGGAGTTATTCGAATCCAAAGCACTATAGAGAGTAGATTTAAAGCCACGCTGAAGTGAAGATAAAATCCGCTTCATCTTTTATGGCGATACACTTTGGTGGCGTATAATGTTGATGCGAGATGAAAAGACCGAACGAATTACTACAGTTTTAAGTTGAGAAAGCATCCACGATGAGCAGTCGGATTTTGGTTGTGGAAGATGAAGCCCGCATTGCCGATTTTCTGGTGCGGGGGTTATCAGAAGAGGGATATGGAGTCGAACACGCGGAAGATGGTCGAATCGGCTGGCTGCTGTTGAATTCGGAGACATGGGACCTGGTGATTCTGGACTGGTGGCTACCGGTGGAAGACGGCATCCAGTTACTCCGTCGGTTTCGACAAAAGAACCGCACGACTCCGGTTTTGTTTTTGACCGCGCGCGATGCGGTTACCGAACGAGTAACTGGACTGGATGCGGGAGCTGACGATTACCTGACCAAGCCATTCGCGTTTGATGAACTGCTGGCCCGGGTGCGGGCGCTGTTACGACGACCCGGCCAGGCGAGCGGCGTGGTGCTGGAATATGCCGACATCCGCGCGGATCTGGAAACACAGCGGGTCACGCGCGGCGGAACGCCGCTAGATCTGACGGCGAAAGAGTTTTCCCTGCTGCTGTTTTTTCTGAGGAATCCGGGCAAGGTACTCTCACGGACGCGGATCTACGATACTGTCTGGGACGATTCGTACGACGGGCTTTCCAATACGCTGGAAGTACATGTCAAGGAACTGCGTCATAAACTGGAAGCCTGCGGCCCCCGCGTGATTCAGACATTAAGAGGACGCGGTTATGTGCTGGAAACGACCTCGGATGATGAGGGAGCCTTATGAAACTGACGACGCGTGTTTCAGCGTTTTTCCTGTCTGCACTCGCGGTGATTCTGATCGGTAATTCTCTGTTGCTATACGGGGTGGCGCGGTCATATCTGAAACATCACTTCGATGAACAGCTCGATTCGCTGCTGCATATACTGGTGGCCGCTGCGGAAGTGGAGATCGACGACGTCAAATTCGAATCCACGGATCATTTTGTGATTCAGGATGCCTATGGAGATCCGGATGATATTTTCTGGCTGGTTCTCTCCGAAGAGGGTCAGGTTGTCGCGCACTCGGATAATTACCACACTGCTTCGGGAAAAACTGTAACAGAGAATGTCAGCTCTCTGGTGAATGCGAGTGAAATCACAAAACCGGGCTGGCGATTAGTTCGACATCATCTGGCCGCGCCCGAACCTAAGCCAGCATCTGAACGTTCGGATTTGGAACACGCCGCGTTAACCATTGTGGTAGCGCGGAACCTGGAACCGTTGCGTCGAGCCCTGTTCTGGCTGGCGGTGGCGTTGATTGTGCTCCCCTTAATCTGCTGGTTAATTGCGGCGCTGCTGGGACGCCGGTTCTGCGAGCGGGCACTGAAACCGATCCGGCAGATGGCAGATGAAGTCCGCATAATTGACGTCCACGACACTCGCGCGCGACTGGACGTGCAGCCTACCCGGGATGAATTGGAAGAGCTGGGTGTGACGTTCAACGAACTGCTCGATCAGCTGTTTCAGGAATACGAACATCAACGGCGATTTGCGGGAAACACAGCGCATCAGTTGAGAACTCCCTTAACTGTGATGCAGGGACAGGTTGATGTGGCACTCAGACGGCCAAGGAGCGCAGAAGAATACCAGGAGACACTGACGACGGTCAGCCAGGCGACTACTTCGTTAAGCCAGACAGTGGAAGCTTTACTCTTCCTGGCACGTCCGGCGGAGGATCAGCCGATTCCCGATTATCAACGGACCGATCTCAGCCTCTGGTTGCCGGAGTATCTGGAACGCTGGAAGAACTCACCGCGGTGGAATGATATCCATCTGAAAACGGCAGCGAGTCTGGTGAGTGAAACCTCACCGACACTGTTGGCGCAGATCATGGAGATACTGATTTCCAATGCGTTGAAATACAGTGAGTCGGGAACGGCTGTAGAAATTCTGGTGCGACGGGAAGCTGCTTTGATTGTTCTGGAAGCCAGAGACCAGGGAATGGGAATTGCAGCGGAAGACCGGGAAGCGATTTTTGAACCATTCTTTCGCACCCGCCAGGCGCGGCAGCAGGCATCGCCGGGTACGGGGCTGGGACTGGCTTTGGCCCGGCACATTGCGACGGCGCTGGAGGGGCGTCTGGTTTGTGTCGATGGTCCGGGAATGGAAACGCGGTTTCAACTTTCACTACCCGCTGAGGGTATCCTTTCCACAAGCTGAAGACTTTTTATCACAGCCTACAGTATCTAATTAAATCACGGCAGCGAGCGCGTTGCCTGAAGTAAGGCTTGCCCCTGCCACCTTTTGCTTCAGGGTCTTCTCTCGTGAGAAGTACGATTTTCTTCTTGCTCCCCATCGCACCACTTCATGTTTTTGTTGAGATGAAGAAAATTTTATCTTCCCTTCAGGGTGCGTTTACGCTGCGCTGCGATACTGCTGATAAAGACTTACATTCCATTGGGATTCGCTGGAGATAATAAAATGACAAAAATTGTATTACTTGTTCAATTAACGGCAGTCGTTACCGCTTTCTTACTAATGTCAACGTGCATCGTTGGAGCTGTAAGGTACAGGTCCGAACAAAAGAAAACTTTGGCGGGGGGCTCATTTCGACCATCTGCTGATAGCCGTTCGACCACAAAGCTGTGTCGTGACATTATTAATTCGCAGTCGTCATAAAGCTCTTCGAATGCAAATATTTCCATGTAGTGTGACGTAGAACGACGACGGTTCTGCGGTATGTGGGACACAGGAACAATCTGTGCAATAGGGAATCGTAGCGATATATATATGGCAGAAATTCATGCTATAAAATGTAAAATAGGAATTCATTTCAAGCTTGGATCGCTTCCAAGTTCAAGCCGTCTAGGAAATGTCTGGTAGTGACACTACTCAAGGCAAAGGGGTCTACTAATTCTCAAAACTCATACAGGCTTCATAGAATCACACTTTACAAACCCTGATAGCCCTAAACTATAATGTTGAATTGACCTTCCTAAAACTCACCAGAAAACCAATAACGATTTCAATCATAGAGGTCCTCATGAAACTCGCTGTTTCTTCCATACTGCTGGTACTGATACTACTGGTCAGTAACTCTCCACTCGCTGCGAAAAAACCACAGGCCGATCATATTCGAGAACTTCAGACCACAGCCATTAAGAATAAAAAGAGCCCGGCCGCACACTGGGGCTTTGATCCGAATAACTATACCCAGTGGTCCAGCCATTCTCTGCGACTGATTCCGGTTTACACTTTCGGAACCCAAAACGGCGTTCACGGTTGCAATTTGGATTCTTATATCGGCAAGAACAGCCCGTATCGCGATGAGAAAAAACTCGAAGCTATTTATGGTTTTCTGCCTGAGAATACGCTGAATCTAAAAGCCAAATATCTGGATCAAACCAATCTCTACGATATCCAGAAGGCAGCCCTCAAAGCAGGTAAAAAGAATATCATCCTGGTTGTGTTCGATGGCATGGATTGGGATACGACCCGAGCCGCGGCCCTGTATTACAATGGGGCTGACAAATACAAAAGTGGACGCGGAACCGGTTTGCATTTTCAGGATTACACGGCTGATGGAACTTCCCAGTTTGGTTTCATGGTGACAGCGCCTCACAACGATGGCTCCAATGTGGACGTCAACACTCAGAAAGTATTGAACCCGGGCGGAAAAATGCGTGGCGGATACAATGCGAAAAAAGGAGGGCCCACACCCTGGAAAGCAGGCGAAGATAAAAAATACCTGATTGGCAGTTCCAGTAATAAATATGGCGAACACGCTTACCCCGATTCTGCTAATACAGCTTCTTCCATGACGACCGGTATCAAATCATACAACAACGCGATCAATGTGGATCCAAACGGCGCTCCTGTCGCGACGATTGCACACGAGGCACAGGATAAAGGATACTCTGTTGGCGTCGTCACCAGCGTTCCGATCACACATGCCACACCTGCCGCCGCTTACGCACATAATGTCAGCCGCAATGACTATCAGGATCTGGCCCGGGACCTGCTGGGACAGACTTCTATTTCACATCCTGAAGAAGCGATGTCAGGTTTGGATGTCGTGCTGGGAGGAGGCTTCGGCACAATGGAAAAACCGACCGGCGGTAAGTCGCATGGCAAGAACTTCGTGCCCGGCTGGAAATACATTTCGGAAGAAACCATCAACAAGGCAGACGTCAAACAGGGTGGAAAATACACGGTCGCTTTGAGAACTCCGCAGGTTAAAGGGAAAGTCGGTCTGAAACAGGCCACTGCCGCCGCCATCAAAAATAAAACGCGACTGCTGGGCGTGTATGGCGTCGAGAAATACGCGGCTCATCTGCCGTTCCAGACCGCTGACGGCGATTACCAGCCTGCTCCGGGGAAAAAAAACAGTGCCGAGGTTTACTCTGAATCAGATATCGCTGAAAACCCGACCCTGGCTGATATGACTGAATCGGCACTCTCTGTACTCAGTCAGAACAAACAGGGTTTCTGGCTTCTGGTGGAATCCGGAGATGTCGACTGGTCTAATCACGATAACAATCTGGACAACTCTATCGGTGCGGTCAAAAGCGGCGATCACGCGTTTAAAGTCATCACAGAATGGGTGGAGAAAAACAGTAACTGGGATGAGACCATCGTGATTCTCACCGCCGATCATGGGCACTATCTGAATATCGACCAGCCGGAAGCGTTCATCCCGCCGAAAAAAGAGGCGAAGTAGAATCCTCTAGATCGCATCACTGATACTCGTAGCGGTAGCCACGATGTGTAATTTAATTCCTGCTTTCTGTGCTGGTGATCTGATCGACGAGTCTGCCGCAGGCTACCGTGCAGTTTATCACGACATTCCCGACATTTAAGCGCGAGTTTTCGCCTCCCACGTCGTAAAATGTCGATGATACCTTCCCACCGTTGTCCTGATTCACAAGAATCAGGCTATGTACAAAAACACGCGCAAGCAACGCAGCTATGATCATCGGCTAAAAGAACTTGTTCGATCTACCGGCAATATCGAATCGTGGTGGCGCGTGCTCAAGCACCAATGGCTTTATTTGAACGAGCTGGATTCCGCACATGCGGTCGAAAACCTGGTGGATTTCTACGTCGAGCAATACAACACTCATTTACCACATTCCGCGTTCCAGGGGCAGACTCCTGACGAGATGTACTATGGAACGGGGCGGGAGATACCGGGGCAACTGCAGGAAACGCGTATCGCAGCGAGAAAATCACGGATGGAGTCAAATCGATCTCAGAGCTGTCGGATCTGTGAGGAACTTCTCGCAGTCGGTAGTTGAGGATGAGAGTGGCAGGACGGGTCTACAGGTTAGTCGTTCGCAGTTTATGATCGACACCGACGATTTCACACTGCGATCGTTGCCTCAGCTTGTAAGTCTTATTTGGTGCGGCGATGATGGTTTGCGTTTCCTTGACTAAGGGGCTGGGCTTTTCCCAGGTGCGAGAACTCGCGCAAAAATGTCTGGAATGTCATGACAAACAGCAGGGTAAGAAATCTTTGGTTTTGCATCCAGAGAAAGAAAATCCATTCACACTCTAAACCGTTTTAGCCCACTCCCTCAAAAATTGCAAAGGAAGATTTCATTCTCGGATAGACTCCTAGACGAAGTTTTAGATGCAAATATAATCCAATCAGGTCCACACAAGGGGAAAGTTCTTACCGCTGATGGTGAGATTATATCACGCTCGGATCCAAGAATTACGATAGAACACATAACGCCTGTTGTTGGACATTGGAATACTGTTGGATACGACTCTACGCGAGCCGTTCGTAAAGAATTTTACAATGACACATCGAACATGTCGATCCGTTTACGTAGTGCTAATAGTTCTGATGGTGGTCGAATGTCTGCAAATGGTATACGCTATCGGCAAGATATTGGCTCAAATTACGAAAGATAGACAGGAGCAAATTTGTGATTTCTCAATCAGACATCGAAGCAAATTTTTTTGATATCGAGTCTATTCCAAGCGATATTTTAAAATTCCAGGATATCGTTCCTGCTTTTAAATTTCTTTATTACTTCATGCCAGACCCTCTTTCACAGAGGGATTTGCTTATAGTGACTGATGGTTATCACGATGATCAATTCGCACAGTCCGTTTTGCTGGCAGTACAACAGTTAGAAGGAGGCCAAAGTGTCTTATCGAAACATCAGATTGGCATTCCAAACGAATACAACTTTACCCATTTGTTATTGGTTCCCTCCGATTTCCACAGCTACTTTAAGGGAAGATTAGATGAGGAACGCAAGGAACTCTATCTTGTTCTTCCGATCCATAATTGTGAATTTTCTGGAAATGAACCGCAGGAATTGTTTGTTCAGATGAGACGACAAACGAACTCCGCTTTGGATTGGCGAAGAGAGCTAACACCGAAGGCCCTGCTCAGATTTGAAAATCCGAGTACACAAGGAGGTGCTGGCAACTCGAATGGTGTCCCCGTACGATTCACACTGATTGATCAGGAAATACGAAATCTTAACGGCATTGAATCAGGTTTTATGGAAGTGACTAGTTTTCGTGATGATTATGTTGAGATACTGTCACCAAAACGGGATGAATACACTTTTCGCAGTCAGTATGACGATGAAGCTCGAACGATGAGTCAAGGAAAAGTGGTCAGTGCCATTTGGGAGTTTCTCGTGGAAAATAATTTGGGATGAGATGTTTTACTTTTTTATTGATTACATACATTCACATATTGGGTCTAGAGCGCGTCACATTTTAGTCTAGCGCGAATGCGGTAAGTCTGGTAGACTGGTGAAGTCTCATTTCATGGAGGATTTCACATGTTGATTTATTCCAAAGAACGCCGCCTCGAAGTTCTGAAAGCGTATGCAGCCGGTTTAACAACCCGGGAGATTGCATTGCAATTTCAATGCAGTGAATCGTGGGTTAGACGGGTCAAACAAGAGTTTCGAGATCAGGGAAAAACCAGCCCTGCCACCAGGCGTAAACGAGTTCCCCAGTGGCATTCACTGGCAGATCGAATTCAAACTGCAGTCTCAAATAAACCGGATATCACTCTGCAGGAATTGAAGGATCAGCTTGGCACAGCACTCTGTCGTCAGACATTGTGTCGCGCATTAACACGTTTAAAGCTCACTCTCAAAAAAAAGTCCTGATCGCATCAGAGCAGGATCGCCCCGATCAGATCACCAGCACAGAAAGCAGGAACTACATACACACTGTGGCTACCGCTACGGTTATCTGTGACGCGCTCTAGGTAAAAGTTCATAGATCGATTCAGTCTATTAGTTCACAATCTTAGCCTGGGAATTTTGAGGACTGGTAGTATGTGTACCCACTATGAATGCATATTAATATTGAGCATCAAAAAACATACCCCAAAAGTCTTCTATACTTGTGATGCAGAAGGTAGAGTTCAGAAGAATAAATTTATCTTGATGTTAAATATATCGAATCGAATTCATGCCTGAAATCACTCGAAAACTAGCGATCTTCAGCTGCCTGCTGTTCACCAACTTCTGTTTCTGGAAAGCGGGTATGGGGCTGACTGGTGGTGTGGATGCCGAGATTCTGCCTTCTGCTCATGCTGCGATGAGCCCTTCCATTTCACCGCCGCGAGTGGTGACGACGCCGATTCAGGAGATGCGGCCCGGGATGCGGGTAGTCGGGCGGAATCCGCTTCGCATTGAAACGGAAGCGACGATCGATCCGACGCCGGCAGGCTGGCGGCTGGTGTCGGTGCGGATGCTCAAACCGGATGGCACTTATTTTGAAGCCGAACTGCTGCGTCCGTTAAGCTGGATTTCGCGGCACCGTGCAAAACCGGGCGCGGTCATCCAGCTCAACATGCCGGAAATGTATGTGGTGGGGGCGGCGGAAGTGCTCTCGATATCAGACTGTCCGCCCATCGATCCGGGTGACGGTCCGGTTGTTATCTCCACGTTCAAAAACACCGCGGACAATGTGCTCAATATCTATGTGGAAGGGGAAACCGAACCGATCGGCGTGACCGCCGGACATCCGATCTGGAGCGAAGACCGCCAGGCATTCATCCATTCCGACCAGTTGCAGCCCGGCGAACGACTGCGTTCTGCGGTGGGTAAGACAGTCCGCATCACGTCAATCGAAATCCGCGCCGGACCGGAACCGGTTTACAACCTGGAAATTGCTGGCGAACACGTCTACAGCGTCACCGGCTCCGGCCTGCTGGTCCATAATGCAGGTCCTTGTTATCTGGTCCCCCGTGGACAAAACGTTGATTTAGACTCGCTCTATTCTACGACAGATGTCCTATCGAGCAGTCATCTTGTTAGTCAATTTGATACATTAAATGCCCGAAGTTTTCACCCAGCGTTAGGAAATCCAAGTCGAGTTTTGACACAGGCTGATTTGAATGTTGCACGACTGCCAAGACGTTACGGCCCTCAAGCAGGCGAAGTTGCGCCATCTAATTGGGGGCAGCACATTATAAGCGAGACTGGCGTTCTTCCGCCCGCTGGTATGCCGAGGTCCCACGGGCATCACATTAACATGAAAGCAGGTCATGGAAGCCAAGTTGAATTTGTCGAAAAAGGGAAAGACATACTTGAGTTTTATGATATCCCTTGGTTTCGAGGCACGGGGCCAACAGGTAACCTCGTTTATGCTCCAAATGTTGCAGGACAGCACACAACAGAAAACGCAACCAAGCTCTACAATGAGCTTCTTGGTGTACATCGCTCAAACATTGACGCTGTACTTACATGGCAACAGGGGCGGGAACTGATTATGGAACAACTTCAAGATGCAGGACGACGCATGTCTCGTCTTGAATTCTGATTCCTTCTGGAGAAGTAATGACAAGTGAAGCTATAATGTATGGAGCCATAGTTCAAGGCGACATAATGGCAGTAAAAGAGTTGGTACGGAACGACCCCTCTATATTGCAAGTGTCAAAAGTTGGAAAGAATTGGCTGCATTGGGCTGCTCAACGAGGACACACCGATATTGCAGCGGTTTTGGTTGAAGCTGGATTGGAAGTCGATAAGTTGACAGACGATGGTACGAGCTCTGCATTGGACATTGCGGCAGGACAGGGAAGATTGGATTCATGCAAGTGGCTCATCGCACAAAGTGCAGAAATCAATCGAGGTTTCGGCAAATGTGCGACTCCAATCTTCAGTGCAATCTATGGTAAATCACTTGAAGTGGTCAAACTGTTTGTCGAAGAAGGCGCACGTTTGGACGCAGAATTTGGCGAACCCGTCATCAATGTAGTTGGTTACGCAAAACGGTACGGCACTCCCGAAATTGTTGAATTCCTTCAACAAAGAACTTCACAACATCCATCAACACCGTCTTTACCTAACGAACCATAGAAGCCCTCAACTCTTCCCTTCTCCCCAAAAAAACGAGCAATTCTAAACCGCTCGCAAGACAATAGAAAATGGTGTCTAAATTACCAAGTGTGAATCCGAAACAGTTAATGGTACCATTCATACCAGAAGTCCTGTTCCGGTACATAGAATTCTGACAAGGTAGAAGACATGAACTTACAAGACATCGTAATACAGGGCGAAGAAATCGTCATTGATGACAGCAAAACATTCAATGCAATCGGAAGTGACGTTGTACTGGAAAATTGCACTATTAGGTGTAGCGTGCCAGCCAAAAGTATGTCAATTCGTGGAAAGTAGGTCAGTTCTCTAGTGATTGCTGAGAGTGATCTTATTGGTTTCTCTTGGCAGGATGCAAACTTGAATCGCTGTACTTTTTAGGGCGTTTTTAGTGATGCTAATTTGGAAGACTGTAACTTTTTCGGTAAATCCTGTGAGTCTCACATGTTCCCCAAGTGGCCAAATGTTGTTATTCTGAATCTGCATAGAGGACATCTCAAAACCATTTTAATAGTGTAAATAAACAACCAAGAATTACAAAGCTTTCGTAGAGGTAATCGTGATATTCCCAGCGCGTGACAATGCGACGATAGTTGTGGAGCCAGGCAATACTGCGTTCCACAATCCAACGTCGTTTGAAGCGTGGGAGCTTTCGACCATCTTGCGTCGGTGGTTTGACTCTCGACTTACGATGCGGACAGATCAGATCGATATTATTTTCCATCAGCCGTTTGCGCAACGGATCCGAATCGGCGGCTTTGTCATAAACCAGTCGCTCGGGTTGTCGATTTTGCAATGTGTTTTTTGCAAGCAGCGGCTCGATCAACTTGACTTCGCTACGGCTGGCCGATTCTGTGTCGATCGCCACAGGTGTCCCGCGACGATCGACGAAGATCATGATTTTAGTGCCTTTGCCACGACGAGTCGGGCCAACTCTTCTACCCCCTTTTTTGCCGAGGCAAACGTGCCATCAGCAAAGGTTTCCGAGAAATCAATTTGGCCGGCATCATCCATTCGTTCTAAGATGATTTGCCAGGCAGACAAGAGCCGTCCTTCGATCGTCCATTGCTGGAAACGTCGATGGCACGTCGCTTTAGAGGGATACTCTTTTGGTAAATCTTTCCATCGCGCTCCTGTCACCAGAATCCAGAGAATGCCTTCCAGGCAATCTCGTGGGTGGGCCTTGGGACGTCCTCCTTTTTTGGAAGGGGGAGTCCAAGGAAACAAATTTGCGACTAAAGACCATTGTTTATCGGTCAAAACGACCGGGCGTTCCGTCCTGGACGCTGTTTTTGTGGAAACTCATTTCGTTTGGGCCACCAGACCAGCGTCATATACATGAGAAATTCTCCTTTCAGGAGAACATTACATGCAAACATCACGCCAAACAGTTCACTTTTTTTGACGTTATGAGACAGCTTCTAGTCACTTTGCCGAAATGCAAGGTCGTCCAAAGCATGTTCGAATTGCAGATGTCGTCGATTCGATTGAGTTTCTTGATGAAGAAGTCAGTGCAATAACTTACAACTCAGAGTCTCTTGCTAAACGACTTGGCCTAAGTATTGACGAAATTCAGAATTTCTTTTCTCAGTTTGCCTTCGTGCAAATGTGAGTCTGATTATTAATTGCTTTGGCTTTTGAATCATGATCATCTCAGGCTGTTCCGTACTGAAAACTTAAATACCAGCTTAACTGCAGATTAAACTCTGCTGGACAATTTATTATGTTATTTGAATCGTATTCATGTCTGAAATCACTCGAAAACTAGCTAGTGTCCGACCCAAAAGCGTTTTGGTTCGTAAATAATAGCGTCAAATTGTAAGAACGCTCCTCGAAGATGAGGCATTGATTCGTATTACTTGTTGAAAAGACCAGTAAAACGGAATCAACACATTCAATCTTCGAGGAGCCAATGATGCGAAAATCATACTCGAATCAGCTGCGACTGGACAGCGTTCCGATCGAACAGGTAGCGTTAAATCTGGAATCGCGGGACCGCATCGTTCCCATTCTTCGCGCACTGCAGTTCCTGTATTCGGATCGCAGGCTTGTCGATGAAATACTGCAGTGGATTGCCGCCGACGTTAACAGCGACAGTCGCACTGACACAGGTCGAACGGGGATGGAATACTGGCATATCTGTGTGCTGGTAGCCGTGCGACTGGGCTGCAACTTCACCTTCGATCAGCTGCAGGATCTTGCTGAAAATCACCGCAAGCTGCGCGGCATCATGGGAGTCGGCGACTGCGATGACAGGCCCTTCAAATGGCGAACCATCCGCAACAACATTCGACTCCTGCGACCGGAAACCATTGCACGTATCAACCAGGCCATCGTCAGTGCAGGACACACGTTTGATCCCACGGCGATCGAAAAAGTGCGGGCCGATTCGTTCGTGATGGAAACGAACATTCATTATCCCACGGAAAGCAGTTTGCTGTACGACGGGCTGCGGAAAATCATTCCTCAGTGTGTAAAACTGGCCGAAGCACACGGCGTGAACGGATGGCGTCAGTATGCTCACCTCTTGAAAAAAATCAAACAACTCAATCGAGACATCAATCGCATCGCGACAAAGAAAGGCCCCCGCTACAAGAAGCGACTCCAGCCGCTTTACCGCGAACTCCTGCAGAAAGCGGCCATCTTGACGCAGCGTGCCCGCGACCTCTGCCTGGTCACCGGCCAGCGGCTGCCGGAGACCGCCGACCTGTTCGGGCCGAACACGCTGCAGGCGTTGATCGTACGCACAGAGCAAGTCGCGGACACGACCCGGCGACGTATCATCCACGGCGAAACAGTTGCCAACAGCGACAAGCTCTTCAGTATTTTCGAGCCGCACACTCAGCTTTACAAACGTGGCAAGGCCGGTCAGCCGATGCAGTTTGGTCGACAGGTTCTGATCTTCGAAGATGCGGCTGGCTTCGTTGTGCGCGCCGTTTTAATGAAACGCAATGAAGGCGACAAACAAGTGGCAGTCCGTGAGACGAAGTCTTTGCAAAATGATTTTCAAAACGGCGTGAAGCGCCTGTCATTCGACCGTGGATTTCATTCTCTCGACAACCAGAGGGAACTCGCCGAGCTTGTCGATCACCTGTGTCTTCCCAGACCCGGCGTCAAGCAGTCGGCGGTTCAACAGACCGATGACGAGTTTCGGTCAGCTCAGCAGAATCATTCGGGAGTGGAATCAATGATCGGAGCATTGCAAAGTGGGAACGCGATGAAGCGATGTCGCGACCGTTCGGAGATCGGTTTTGAACGCTACCTGCAACTGGGCATTCTGGGGCGAAACCTGCACACGCTCGGTCGGATGCTCATCGCCCGGGAGAATCAAAACGCGGCCGCCGCTCACAGTCGCCGCAAAGCGGCCTGAGACAGTCGCGGTCGATACGAGAGAGCAAAACCAGGGAAGAACGTACGCGCGGCATTGTCTGATGAAGTCGAAATCCCCGCCAAAATCACCCGAACGCAGCCAGCAGAGACCCGATAGCGACAATCTCCACGATCCACTACCGCCGCGTTACCACGTAAAATTCTGTACGAAGCTACTTTCGGGACACCCACTAGCTATCTTCAGTTGTTTCCTGTTCACTGGTTTCTGTTTCTGGAAAGCGGGAACGGGACTGATTGGTGGTTCGGATACCGAGATTCTGCCTTCTGCGCAGGCAGCAATGAGTCCCTCAGTTAAACCGCCGCGGGTGGTGACGACGCCGATTCAGCAGATGCGACCGGGTATGCGGATTGTCGGTCGGAATCCGCTTCGCATTCAAACAGAACCGACAATCGATCCGACGCCGGAAGGGTGGCGGCTGGTTTCGGTGCGGATGCTGAAAGAGGATGGGACTTACTTCGAAGCAGAACTGCTGCGTCCCTTAAGCTGGATTTATCGACACCACGCGCAACCAGGGGCGGTCATTCAGCTCAACATGCCGGAAATGTATGTGGTCGGGGCGGCGGAAGTGCTCTCGATATCGGACTGTCCGCCCATCGATCCGGGTGACGGTCCGGTGGTCATCTCCACATTTAAAAACATTTCGCATAACGTCCTCAACATTTATGTGGAAGGGGAAACCGAGCCGATCGGCGTGACCGCCGGACATCCGATCTGGAGCGAAGACCGCCAGGCATTCATCCATTCCGACCAGCTGCAGCCCGGCGAACGTTTACGTTCGGCGGTGGGTAAG";//ACAGTCCGCATCACCTCCATCGAAATCCGCGCCGGACCGGAACCGGTTTACAACCTGGAAATTGCCGGCGAACACGTCTACAGCGTCACCGACTCCGGCCTGCTGGTCCATAATGCAGGTCCCTGTTATCTGGTCCCCCGCGGACAAAACGTTGATTTAAATGCTCTAAAGGTTGCAGACGGACCATATGATTCGCGCACAATTAGGTCGCAGCTTGAAGCTCGATATGGAGCAGAGGCAGTCACATCTACAACAGTGCCACCCACCTCTGGAAAAATGGTCCAAATGTCAGGTCGTCAACACTCTGTCACGGGAGTTGTTTATGATTCAAAGGGGTTTCCTATCTTTGATCGTCATTCATTATTCGACACATATATTACACGACAATACAGCAGTGTGGCAGACGAGGCAGCACATATGCGAGCAGCAACCCGTGAACTCCGATCTAGTTTGGGCATCGGACAGGGTCAAGCTAATAGAGTCGCTGCGCTCGGAGCGATTAAGAATAAAGTTAACTTTCAAGGTTTCAATCCTGTACAGCTTCGTCAAATTATGGGGGGAAAGAAAAATATTGATGGTTTAACTTGGCACCACCATCAAGATATTGGGAGAATGCAACTTGTAGACCGCGTGATTCATGCGGAAACTGGGCATGTTGGTGGTTTCAAGATGTGGTT";//TGGATACTAATAGGCAAAAAATGGTTAGTATAAAAGAAATTGAATGGAGCGAATGCAGAGGCAAACCTACGGTGGATGCGATAAAACAATTCGAAGTTGTGATAGGTCGCGAATTGCCTGAAACTCTTAAGGTACTCCTAGATACAGTAAATGGTGGAAGTCCGTCTCATGACATTTTTTCTTATAAAGATCCCGAAACCTATCAAGAGATCCGTTCTTGTCTGGGAGCGTTGATCAGTTTCAATTCAGAGGATGAAAACAACATCTTAGATTATCTCAGCCGACCTTCCGAAGGAATGTCTAGTGAATTAGCGCCATTTGCAGCGGTCGGAAACGGTAACTTGATCTGTTTAAAAAAAAATAGCGAATGTGTCTTCATCTGGCTGCATGGAAATAGATCTGGATATGACGAGTGTCAATTAGCTGATACATTTCAGGAGTTTATAGAAGTGCTTCAGCAAGATAATGATGATGAAGACGAAGACGA";//AGACTTGGATACCGAGTTAGCAAATCTAGAGGAAGAATAGATAGATAAAAAACGTACTATTCAAGACAAGAAGCTTAAAATATTAGAAAATAATTCACAAGGATAATTAAAAAGAGTTAGCTGATGAGGGGGGGTAATCTTTGATGAATAAAGACATTCTTTCTATCACAATCACTCATTAGTTTTTTGAAAACAGGAAAATTCCCACAATCATGACTGCTCCATAAAACAAAAAAATTATTCATGTCTGACCTCACTCGAAAACTAGCTATTTTCAGCTGCCTGCTGTTCACCGGTTTCTGCTTCTGGAAAGCGGGAATGGGTCTGATTGGTGGTGCGGATACTGAGATTCTGCCCTCGGCGCAGGCGGCAGTTCACCCATCCTACATACCGCCGCGGGTGGTGACGACGCCGATTCAGCAGATGCGACCGGGGATGCGGATTGTCGGTCGGAATCCGCTTCGCATTCAAACAGAACCGATAATCGATCCGACGCCGGAAGGGTGGCGGCTGGTTTCGGTGCGGATGCTGAAAGAGGATGGCACTTACTTCGAAGCCGAACTGCTGCGGCCCTTAAGCTGGATTTATCGACACCACGCGCAGCCGGGCGCGGTCATCCAGCTCAACATTCCGGAACTGTATGTGGTGGGGGCGGCGGAAGTGCTCTCGATATCGGACTGTCCGCCCATCGATCCGGGTGACGGGCCGGTGGTGCTTTCCACATTTAAGAACATCGCGGACAATATCCTCAACATTTATGTGGAGGGGGAAACCGAGCCGATCGGCGTGACCGCCGGGCATCCGATCTGGAGTGAAGATCGGCAGGCGTTTATCCAATCCGACCAACTGCAGCCCGGAGAGCGCCTGCGTTCGGCGCTTGGTAAGACAGTCCGCATCACCTCGATCGAAATCCGCGCCGGACCCGAACCAGTTTACAACCTGGAAGTCGCCGGCGAACACGTCTACAGCGTCACCGGTTCCGGTTTGCTGGTCCATAACGCAGGTCCTTGTCGCCTGACTTCAAGGGGGCATTACGATTCGGGCTATATCGGCATCGTCAATGGTTCCGCTGATGACTTCTCAGATACGCTCCGCGCCACCTTCGGAACGTCTGACGTGCTGCCAAGCAGTCGACTTGTCTCTCAGGCCGATACGGTTGGCGGGATTGGCTTTAGCGCCTCAAGGGGGGACTCTATGATGACATTTCATTCCCAAGCGGTCGCTACAATGATGCAACACTGGATAAGTATGTTCCAGGTGGAGCCGACTATGCCATTTTTCGAATGGATGAGCGCTTGATGGTTAATATGCCAAGAGCCTTACCTGGTCAGAGCCAAAATGACGCAGGATACCTGCGAGACGCTCAATACTACTGGAAGGAACTCTACAAAAAATACCCAGAAGCATTTAGTGAGCAAAATATCAAAATAATGAATGGACTTGTGAAGAGCAGAATTGGTAATGTAGTCACAGCGATAAAAAATGACGCAACGTTTAGAGCTGTTTTTAAACAATACGATGTCAAATACCTAAGAGGTAAATCAATGGTTCACCATCATGTGGGGGGGAAGCATAGCTGGTGCAATACCATCACCTCTTCACCCTGGATCTGGAGGAATCCATAACGTTGAAAAGCAACTTGGAATTATTTATGACTAGGTAGATCGAATGACCACAAACAATACTTTCGTCGCATTTGCACGAGATGGCAATCTGGAAGGGATCAAACACCTTGTTTCAGAAGGTGTGGACCCAGACACTACGAATGTCGCAGGTGCAAGTGCCATTGAATTGGCATCACAGAATGGACACGGAGAAATAGTGTCGTGTTTGCTTTCGCAAGGCGCATCTCCTGATCCAACGACTACCACGAGACCTCATTCGTCTGCATTGATTTATGCGTCAACTAATGGTGATCTGGAAATCGTTCAACAACTGATTTTACATGGTGCAAACGTCAATCTGATCGATGATGCCTACGAGCAGACCGCATTGATCTGGGCATGTGCCTATGGCAAGTCTGCTGAAGTAGTAAAAACCCTTCTCGATGCCGGTGCAGATGCATCAATTAAGGATGTCAACGGAAACACAGCTCTAGTACTTTCTGAGCATTTCGGTTTTGACGAGATTACCAGTATTTTATTGCAGAACTCATCTCGATAATCGTGCACCTGCTTCTAGTCATCCAACGAGGAACAGTTAAAGGTGCTATTCAGAACTGCCCCCTTTAATCTCCTCCCAAGTTCGTGAAATCTCGTTGACAGCTTTTGAAATCAACTGCTGAGGGTTGGCAAGCAAGTGAATTCGCTGATTACTGTGATTTACTTCATGCGCAAGTTTACAAACCACTTGCCGATGCGCTTATCTATGGGCTTAAAAAAACTCGGAAGTAGAATTGATAAATCAACTCATCCTCCTAAAAGTTGATCCTCGTTGATGATTTCTCTAAGCGTTGATTCAGTCGTAACACGGTCATGCAGGTTCTGGATCCAAAATACGAAATGTTAGTCAAACTGAGTAAATAAAGTGTAAAAATGGATTTTTGACAAGCAGACAAAATGTTTTTTGTAATTGGATTCGAATTCATGCTTGAAATCACTCGAAAAATAGCGATCTGCGGCTGTCTGCTGTTCACCGGTTTCTGCTTCTGGAAAGCGGGAATGGGTCTGATTGGTGGTGCGGATACTGAAATTCTGCCTTCGGCCCAGGCGGCAATGAGTCCGTCTCTCTCACCGTCGCGGGTGGTGACGCCGCCGATTCAGCAGATGCGACCGGGGATGCGGATTGTCGGTCGGAATCCGCTTCGCATTCAAACAGAACCGATAATCGATCCGACGCCGGAAGGGTGGCGGCTGGTTTCGGTGCGGATGCTGAAAGAGGATGGCACTTACTTCGAAGCCGAACTGCTGCGTCCCTTAAGCTGGATTTATCGACACTACGCGCAGCCAGGGGCGGTCATTCAGCTCAACATGCCGGAAATGTATGTGGTGGGGACGGCGGAAGTGCTCTCGATATCGGACTGTCCGCCCATCGATCCGGGTGACGGGCCGGTGGTCATCTCCACATTTAAAAACGTTTCGCACAACGTCCTTAATCTTTATGTGGAAGGGGAAACCGAGCCGATCGGCGTGACCGCCGGACATCCGATCTGGAGCGAAGACCGCCAGGCGTTTGTGCATACGGACCAGTTACAGCCCGGCGAACATTTGCGTTCGGCGGTGGGTAAGACAGTCCGCATCACCTCCATCGAAATCCGCGCCGGACCGGAACCGGTTTACAACCTGGAAATTGCCGGCGAACACGTCTACAGCGTCACCGATTCCGGCCTGCTGGTCCATAACTCAGGCCCCTGTGACCTCGTCCCTCGCGGGTATTACGATTCGGGTTATATCGGCATCGTCAATGGCTCCGCTGATGAGTTCTCAGATACGCTCCGCGCCACCTTCGGAACGTATGACGTGCTGCCTAGCAGTCGACTTGTCTCTCAGGCCGATACGGTTGGTGGCGTTGGATTGGGTTCCCAAAAGATTGGAACTAGTCAGGGAGTATTCACAAACCGGGTAGAATGGTTGGCTCCAAAGGGGACGGAGCAATTGTATACAGTATTTCAAAGGAACGACATTGACTGGGCTTTCACCCGATCAACCGACTCGGGAGTCCTTGCTGGAAGAGGGCTTACAAACGCCGAAGCGGCTGCAAAACAAGGACTTGGACCAATTCTACCAGATGGTTACATTGCAACTTTGCATCATTCACAGCAAAAGTCTATTGGTCCGCTGTTCAAAGCATCAACTCGCTATCATAACTTTTATAGAATAGATAAGGCACCATTGCATCCAAACGGTTATAGACAGCACCCAGATTTTCCAATGGGTAGTGGTCGTGATCCGAATTCTATCCGTTCTGCTTTCCAGAACCATGATTCGTTGGATTACTGGAAATGGCGGGGATAGCAAGCGTTAGAAGAACTTGGCAATAGGGGTGGATAATGGATAGTCGAATAGAACTCATTGTCGACAACTTGATTTGTTATTGTAACCAGCGTGGTTTGTCGGATAGTTGCCGATCTGGCAACGACCCGACAATAAAGGTGGTGCTAAAGAAGCTTGGAATAAAGGCCGACTCGGCTCTTGCGTATTTGTATCAACATATAGTTGGGCCATTTTCGAATAATGGTGAAAGACAGCTTCCAGAGTTACTCGACGTTAAGAGTGGTGTGCGGAACATTGCGTCATTTTCAGTAGAACTTTGGGAGAACTATTCGCTTCCAAAAACTATTTTGGCGATTAGCGAAATAAATACCATGACCGGTCTATTCTACGATACCGAGAGCGACAGTGTATACCAAGTAGAAATGGACCTTGAATTTCCTTATTTTCTGAATGGAGAAATCTCACCTCAATGGTCGAATTCGGCTGAATTCTTGGTAGATTTTTTCGGAGAGGAGTCAATTGCGAAAATCTAAAAATGCTTTTGATGAGCTCTTATCAATATTACTAACCACATCATACCTGTAAGTCTGGTTCAATCCTACTTGTACCCGAGATACGAGGACGCAAAACGGGGGCCGGTCTCAATTCTGGTGCGTGGTTTCAGGAGAATCAGAAAAAATTCGTTCTTCAGATATCTGGTGCAGAAA";
//		String b = "AT..T.GT.CA.TGG.T..G....AT........G.AT..G.AA..G.CC..G.TT..G.....G.AA..G..G.AA..T..G..G.AC.TGGGC.CA.AT.AA....AA.AA.AT....AC.GC.AA.GT..T.AA.GA..G.GT..T..T.AC.GT.GT.GC.AA....GA.AC.GA.AC.GT.AC.GC.GG..T....AC.GT.AT.CA.GT.AC..T.AA.GA.GT.GA.CC.AA.GC.GG.AA..T.AT..T..T..........AC.CC.AC.GC.CA.GG.AC.AC.....T.GT..G.AA.AT.CC.CC.GC..T.CA.GC.AA.TA.AA.CA.TT..T.AA.TA.ATGCC.GG....AC.AT.GA.GA.TA.AT.TGGATGAA.TA.TGGAA.GG.GG.TGGGA.CA.TA.GC.CA.CA.AC.GG.AA.GA.GT....TT..G.GT.CA.TA.AC....AC.GA.GC.AA.TA.AC.GG.AC.GA....TT.AC.TA.GT.GT.CA.GA.GC.GA..T.AA.CA.GA.AC.GG.AC.GT.AC..T.AC.GT.AA.GC.AA.CA.CC.CC.GT.AC....CC.GA....GC.AC.GT.AC.GC.GG.GA.TT.GT.AA.ATGCA.GT..G.GC.AA.GA.TA.GA.CC.AA.GG.AC.AA..T.GT..T.GA....TT.AC.CA.CC.ATGCA.GG.AC.GT.....T.AA..G.AA.TA.CC.GC.....T..T.GC.GA.TT.AA.CA.TT....AA.GA.ATGCC.GG....AC.AT.GA.GA.TA.AT.TGGATGAA.TA.TGGAA.GG.GG.TGGGA.CA.TA.GC.CA.AA.AC.GG.AA.GA.GT....TT..G.GT.AA.TA.AC....AC.GA.GC.AA.TA.AA.GG.AC.GA.AC.TT.AC.TA.GT.GT.AA.GA.GC.CA.GG.AA.GC....AC.AA.AC.GT.AC.AT.AC.GT.AA.GC.AA..T.CC.CC.GT.GC.AC.GC.GA.AA.GT.AC..T.AT.GC.GG.TT....AC.GT.AT.GA.GT..G.GA.AA.GA....GA.CC.AA.GG.AC.AA..T.AC..T..T....GC....GG.GC.AC.CA.GG.AC.GT.GC..T.AA..G..G.AT.CC.AC.AA..T.CA....AA.TA.GC....AA.GC.TA.ATGGG.CA....AT.AA.GA.TA.GT.TA.GG.CA.TA.GG....TGGAA.AA.TA.GC.TT.TA.ATGCA....AA.TA.TT.AT.AC.GT.GA.AT.GA.AA.TT.GA.GT..G.TA.AC....AA.GA.GC....TA.GT.GG.AC.GA.AC.TT.AC.TA.GT.AT.AA.GA.AA.CA.GG.GT.CA.CA.AC.GG.AC.GT.AC.GT.AC.GT.AA.GC.AA.CA.CC.CC.GT.GC.GT.GA.GA.AC.GC.AA.GT..T.GC.GG.AA.AC.GT..T.AT.AA.GT..G....AA.GA.AC.GA.CC.AA.GG.AC.AA..T.AC..T.AT....GC.GG.CC.GC.GC.CA.GG.AC.AT.GT.GT.AA..G.AA.AT.CC.CC.AA..T.CA.GC.GA.TA.GA.GG....TGGTA.GC....AA.GG.TA.AC.AT.GA.CA.TA.GT.CA.GA.ATGTA.ATGAA.GG.GG.TGGAA.GC.TA.GC.CA.CA.AC.GG.AA.GA.GT.GG.TA.AT.GT.AA.TA.AC....AC.GA.......TA....GG.GT.GA....TT....TA.AT.GT.AA.GA.AA.CA.GG.GT.CA.CA.AC.GG....GT.AC.GT.AA.GT.AT.GG.CC.AC..G.GT..T....GT....GG.GC.AC.AC.GT.CA.AA.GA.CA....GT.AC.TT.AC..T.AC..T.GA....AA.GC.AC.AA.AA.AT.AC.GT..T.AA....AC..T.GC.GG.AC.GC.GG....GA.GA.TA.AC.AA.AT....AA.CA..T.AT.AC.TT.GC.GT.GG.AA.AC.AC.AA.AC.GT.AC.GT.AC.GT.AA.GA.GA.GC..T.GT.AA.GC.GA.AA.GT.TT....TT.GT..T....GA.CC..G....GG.GG.AT....GA.AT....AA.GT.AA.TT....GA.GG.AC.AT.GA.GA..T.AA.AC.CA.GC.AC.AT.AC.GA.AA.GA.AC.GC....AT.AC.GT.AA.GA.GT.AC.GT....AA....GC.GG.GT.GC.AC.TT.GT.GT.AC..T.AA.CA.GA.GT.GA.AT.CC.AT.AC..T....TA.GT.GC.AC.GG.AA.AC.GC.AC.GG.AC.CC.CC.GG.TA.GC.GT.AA.CC.CC.AA.GG.TA.GA.TT.TGG...GT..T.GG.CC.GG.CC.AA.TT..T.CC.CC.AA.TT.AC..T.TA.GG.GC.GG.GA....TT....TT.GA.GT.AA..T.TA.GA.GA.AA..T.GT.AA.GC.GA.AA.GA.TT.AC.AT.AA.TT.GC.TT.GC.AC.GA.GA.TA.AA..G.GA.AT.AT..T.GA.GG.AA.GC.CA.GC.GT.AT.CA.GA.GA.GA....GC.GT..T....GT.GG.GA.GT....AT.GC.AA.GG.GA....GG.GT.AC.....T.AA.TT..T.GT.....T....AA.CC.GT.GA.GA.GC.GT.CA.TT.AC.GC.GT.AC.GC.GA.GG.AC.GC.AC.GT.GC.GA.GG.GA.TA.AC.GC.TT.AA.CA.AA.TA.AT.AT....GC.GG.GA.AA.AA.GT.ATGGT.AC.GT.GA.GT..T.GG.AA.AC.GT....AA..T.GA.AA.AC..T.CA..T.GT..T....AA..T.AA.GC.GG.GG.TT.CC.GT.AC.TT.AA.GG.GG.GG.......AT....GC.AC.GG.AC.AT..T.GA.GA.GA.AT.GC.GC.AC.ATG....T....GA.GC....GT.GT.AA.GG.AA....GG....AC..T..T.AC.TT.AC.GC.AC..T....AA.AC.AT.AA.AA.GA.GT.AC.TT.AA.GC.CA.GT.AT.GA.GG.AC.GC.....G....GA.AA.GA.TT.AC.TT.AT.GA.GG.TT.TA.AC.AT.AT.GC.GG.AA..T.....G.GA.AT.AC.GT.AA.AT....GG.GA.AC.AA.GT.AA....GC.AA.GA..T.CA..T.AT..T....AA..T.ATGGC.GG.GA..T.AA.GT.GT.TT.AA....AA.AC....AC.AA.AC.GC.AC.GG.AC.AT.CA.AA.GA.GA....GC.GT.AT....GT....GA.AT.AA.GC.GC.AA.GG.GT.GC.GG.AA.AC.AC.TT.AC.TT.AA.GC.AC..T.AC.CA.CC.GT.GA.GT.GA..T.AC.TT.GA.GC.GA..T.TT.AA.GG....GC.AC..T....GA.AA.GA.TA.AA.GC.TA....CA.AA.AT.AC.GT..T.GC.GG.AA..T.AC.CA.CA.AT.AC.GT....GT.CC.GG.GA.GT.AA.GT.AA.GC.GA.AA.AC.TT.CA..T.GT..T....AA..T.ATGGC....GG..T.AA.GT....TT.GT.GG.GG.......AC.AA.AC.GC.AC.GC.AC.AT..T.AA.GA.GA.......GT..T....AT.AA.GA.ATGCA.GT.GT.AA.GG.GA....GG.AC.AA....AT.AC.TT.AC.GT.AC..T....GC.CC.GT.GA..G.AA.AT.AC.TT.GA.GT.AC.AC.GT.AA.GG.AC.GC..T.GC.GG.AA.GA.TT.AA.CA.CA.GA.AA.AC..T.AC.AT.AT.GC.GG.AA.AC....AC.CA..T.AC.GT.AC.AT.CC.GG.GA....AA.AT.AA..T.GA.AA.AC.TT.GG.AT.AA..T....AA..T.AA.GC.GG.GA.TT.AA.GT....TT.GC.GG.GG....GC.AC.AA.AC.GC.AC.GG.AC.AT..T.AA.AA.GA.TA.CC.....T.CA..T....AA..T.AA..T.AC.AA....GG.....T....GG....AT.AC.GG.AA.TT.AC.CC.AC.AA.TA.TA.GT.AA..T.GA.TA.GG.CA.AA....GT.GC.GA.AA.AT.GT.GT....AA.AC.AC.GA.TT.CA....CA.AT.AA..T....AA..T.GG..T.GG..G.AC.ATGCA.AT..G.GT..G.GG.AA.AA..T.GC.AC.GG.AA.ATGGT.TA....GA.TGGGT.AA.AC.GA.TT.GT.GC.AT.CC.CA.GT..G.GA..T.AA..T.AT.AA.GA.AC.GG.GT....AA.GA.AA.GT.AC.GC.GA.GG..G.GT..G.GG.GT.GT.AA.TA.CA.GG....AT..G.GC.GA.AA.GT.TA.GT.GC.GT.GA..T.GA.AA.AA.GG.AA.AT.GA.GC.GT.GT.AT.CC.GA..T.AC.GG.GC.TT.GT.AT.GA.CC.GC....TGGGG.GC.GA..T.AC.GT..T....AC.AT.AA.TT.CA.GC.TA....CA.AC.AC.GG..T..G.GG.GA.TGGGT....AT....AT.AA.GT.GT.GA.CC.GC.GT.AT.AA.CC.AA.GG.GC.CC.TGGGC....AT.AA.GA.GT.AA.GA.AT.CA.GT.AA.TT.GG.AA.GT.GT.AA.CA..G.AT.AA.TT.AA.CA..T.GT.AC.GA.CC..T.GG.AA.TA.AC.TT..G.GT.AA.GG.CA.CC.GC.AA.GG.GG.GC.GC.ATGAC.GG.CA.CA.GG.ATGAT.ATGGT.AA.CA.GG.AA..T.TT..G.TGGGA.CC.CA.....T.GA.TA.GG.AA.GC.GC.AA.TA.GC.GG.AA.AC.AC.AT.GC.AA.AT.GG..T.TA..G.GA.GG.GT..T.AT.GA.AC.CA.AC.AT.AA.AT.AC.GT.GG.GA.AT.AA.CC..T.AA..T.AA.GC.AT....GA.TT....GT.GA.GT.TA.GA.CA..T.GC.TT.AC.GC.AA.TA.AA....CA.TA.GA.CC.GC..T.AC.GG.CC.GC.GT.TT.GT.TT....GA.GG.AA.GC.GA.AC.CC.GA.AT.GC..G....AT..T..T.AA..G.GG.GA.GC....GT.AC.TT..G..T.CA.AA.AA.AA..T.TA.TA.TT.GG.CA.AA.GA.GC....CA.AT.TA.GA.TA.AC.GT.AC.AC.TA.TT.CC.CA.GC.GG.GA....GA.CA..T.GG.GC.AA.GC.TT..T....GG.GT.CC.GG.AC.AC.CC.GT.TT.CC.GT.AC....AA.GC.CA..T.AA....AA.GT.............CA....TT.AA.GT....GT.GG.GA.TA.AT.GC.CC.CA.GT.AC.GC.AC.CC.AT.GC.GT.CA.AC.AT..T.AA.GG.AA.AC.TT.AA.TA.GA.GT.AC.GG.CA.TT....CA.CA.TGGGG....GC.CA.AT.TT.CA.AA.AA.GA.GT.GT.AC.TA.AA..T....CA.GC.CC.AA.GG.ATGCA.AT.AA.AA.GG.GT.AT.AC.TGGAA.AC.GA.AT.GC.AA.AC.CC.CC.GG.AA.CA.AT.TT.AA.GT.AT.TA.AA.GT.AA.GG.TA.CA..T.TT.AC.....T.AC.CC.GT.GT.TA.GT.CA.AA.AT.AA.CA....CC.GA.TT.AA.TT.AA.GG.GG.GC.GA.CA.CC.AA.AA.AT.GT.TT.AC.GG.CC.AC.AC.TA.CA.AA.GA.GA.TGGAT.AC.AA.AT..T.GC.......AT.GA.CC.GA.CA.AA.AC.AA.CA.GT..T.AC.TT.GA.AT.CA.AT....GA.GA....GA.GA.AT.TT.GC.GT.GG.GG.AC.CC..G.AT.GT.GA.AC.GT.GA....GG....TA.AC.TT.AC.TA.GA.CC.GT.AC.GG.GG.TA.CA.AA.GG.GC..T.GC.TT.AC..T.AA.GC.GA.GT.GC.GG.AC.GC.AA.AT.ATGAT.AC..T.AA.GA.GG.GG.GG.AC.AA.TA.GG.GG.AT.GA.AA....CA....TA.AC.TT.CA.AT....GC.TA.GA.CC.CC.AA....AC.AA..G.TA.GT.AA.CA.GA.AC.AC.GG.GT..T.AT.....T..T.AA.AA.GA.AT.GC....GG.GA..T....ATG...TT.AC..T.CC.GC.GC....AC.GG.CA.GG.GG.AC.GT....ATGAA.CC.GG.AA.GC.GG.AC.GT..G.TA.AT.CC.CC.GC.GG.TT.AT.GG.CA.GA.AC.TT.AC.TA..T.GT.AC.TA.CC.GA.GG.AC.AA.AA.CC.AC.AA.ATGGT.GT.AA.GT....CC.AC.AA.GA....AT.GC.AC.GC..G.AA.AT.AA..T.GC.TT.AA.AA.AC.GT.GA..T.AA.GG.GG.AT.....G.TT.GA..........GA.GT.GA.TA..T.GG.TT.GA..T.AA.GC.AA.AA.GT.AT.AT.GT.GA.TT.CA.GC.GA.GA.CC....GA.AA.TT..T.GT.....T.GT.GA.AA.AA.GG.GA.CC..T.AA..T.AT.AT.AA.CC.TT.GA.GG.AC.AA.TA.TA.AA.GC.AC.GC.GA.GG....TA.TA..T.AA.AT....AA.GG.AA..T.......AA.ATGAA....TA.AC..T.TT..T.AC..G.CA.GC.AA..T.GC..G.GA.GA.....T.GT.GT.TT.GA.GA..T..T.AA.CA.AT.GA.GT.GC.GA..T..T.GT.GA.GA.GG.TA.AA.AA.GG.GC.GA.GT.AC..T.GT.....T....GC..T.AC.AA.AA....GT.GC.GT.GC.CA.AA.GG.GC.AA.AT.GG.AC.GC.AC.CA.AT....TA....CC.AA.TT.GG..T.AC.GC.CC.CC.GC.GA.AT.CA.GA.GA.GC.TT.TT.TA.CA.GT.CA.GA.GC....GG..G.AT.TA.AC.GG..T.GT.....T.AC.GT.AA.....G.AT.CC.CA.GT.CA.GA..T.CA..T.GT.AT.GA.AC.GT....GA.AC.CA.GG.CC.GG.AA.AA.AC.GA.GG.GT.AC....TA.CC.GT.AT.AC.GG..G.AT..T.GA.CA.GA.CA.............TA.AA.GT.CA.AT.GA.AC.GA.CA.GA.GG.AA.GC.GA.GT.....T....GC.AA.AA.AT.TT....TA.GA..T..G.AC.AA.GA.AT.AT..T.GA.GG.AA.CA.AA.AC.TT.AA.GT..G.GC.GT.AA.GG.AC..T.AA.GG.AA.TGGCA.AC.AT.GT.GT.AC.TA.CA.AC.AC.GG.CA.CC.AA.AT.CC..G.GT.AC.GT.AC.AA.AA.AC..T.CC.GG.GA.GG.AC.AA.AA.GG.GA.AA.GA.TT.AC.TT.AA.AT.GC..T.AA.GA.GG.GC.GT.GC.CC.GC.GG.GG.AT....GT.AA.GC.TA.AT.GA..T.GG.GA.GC.AA.GT.AA.GA.TT....AA.TT.AA.GG.AA.CC..T.AC.AT.AC.GA.GG.AA.GT.GA.CC.GT.AC.GT.AC.AT..T.AA.GG.CA....AA.GC.AA.TT.AC.GT.AC..T..T.GA.GA.TG..T.GT.AA.AC.GC.AA..G.TT.AC.AT.TA.GT....GA.CC.GC.TA.CC.GC.GA..T.AA.GC.AA.AC.GG.AT..T.GG..G.GC.AA.GG.AT.AT.GC.CC.GA.CC.GA.GA....GC.....T.GA.GT....GC.AT....CA.AA.AA.GG.GA.......ATGGG.AA.TT.AC.GA.TT.AT.TT.AC.GT.GT..T.AC.GG.AA.GT.AA.GG.GG.TT.AA..T.CC.TT.AC..T.CA.GA.GG.AC.GC.AC..T.GC.GA.GA.GA.TA..T..T.CC.GA.......AC..T....TT.GC.GG.AC.GC.GG.AA....CA.GA.AT.GT.GT.CA.GT.AA.GC.GA.GA.GT.GT.AA.TGGGA.AA.AA.TT.AT..T.CA..T.GG.AA..T..T.AA..T.GA.GC.....T..T....AA.GT.AC..T.GC....AC....GT.GC.GG.AC.AT.AT.AA.GA.GA.AC.GC.AA..T....AT..T.GA.GA.GT.AC..G.TA.AA.GG.CA....GG.AC.AC.AA.TA.AA.TT.AC.GT.....T.GA.CA.CC.AT.GA.AT.CC.AT....GT.AC.GT....AC.GT.GA.GG.AC.GC.AC.GT....GA.AA.GA.TA.GT.GC.AA....AA.AT.TT.AA.TT.GC.GC.AA.AA.AC.GG.AC..G.AC.TT.AC.GT.AA.GT.AA....GA.AC.AA.TA.AA.AA.GG.GC.GA.AA.AC.TT.CA.GT.GT..T....AA..T.GT.GC.GA.GC..G..T.GA.AA.AT.AC..T.AA.AA.GG.AC.GG.AC.AT..T.AA.GA.GA.GA.CC.GT.TA.AA.GT.AT.AC.GT.AA.GA.TT.AA.GA.GC.GA.AA.AA.AA..T.CC.AC.AT.AC.TT....AA.GC..T.AA.ATGGC....AT.AA.......AA.AT.GA.AT.AT..G.TT.GA.CC.GG.AT.AT....AA.AA.AA.GT.CA.AT.AT..T.AC.AC.CC.AT.AC.AT.GA....AA.GT.AT.AT.GA.GG....GT.CA..G.AA.AA.GA..G..T.GT....GG.AC.GC.AC.CA.AC.GC.GG.GC.TT.GA.GT.GG.GA.TA.AT.CC.GA.GG.AC.AA.AC.TT.AC.GT.AA.....G.GT..G.TT.GT.AT....CC.AA.GC.AC.GA.TT.AA.GG.CC.ATGAT.ATGGT.AA.GC.AC.CA.GG.GA.GA.GA.AA.AC.AA.AC.AC.GC.AA..T.AA.AA..T.AT.AT.GA.GG.AA.AA....GG.GT..T.CC..T....AA.GC.GG.GA..T.AC..T.AT.AA.GT....GT.AC.GC....GT.GA.AA.AC.AA.CC.GG.TT.GT..G.AA.CA.AT.AA.GG....GG.GG.GG.AT.GT.AA.GC.AA.AA.GC.AC..T.AC..T.AA.AA....GT.GT....GG.CA....AT.AC.GG.AA.GG.GG.GC.AT.TA.AA.AC.GC.AC.GC.GT..T.AT..T.AA.......AC.AT.AA.AA.AA.AT.GC.GC.GG.AA.GG.GC.GG.GT.TT.AA.AA....GG.GG.AA.AT.CA.AT.AA..G....AC.GT....GG.AA.AT....AA.AC.GG.AA.GG.GG.GG.AT.TA.AA.GC.GG.GT..T.AT....AA.......AC.AT....GG.AA.AC.GC.GA.AC.AA.GG.GG.GG.AT.TA.AA.GC.GC....GG.AA..T.AC.GT.AC....AA.AC.AT.GC.TT.AA.AC.GT.AC.AC.GG.AC.GG.GC.GG.AT.TA.AA.GG.AA.GT.AC.AA.GA.AC.GC..G..T.CA.AA.AC.AT.GT.GC.AA.AA....GA.AA.CA.GA.AT....GG.AA.TT.AA....AA.GG.TA.AA..T.AT.GG.AT.GT.GA.AA.AA.CC.GC.GG.AA....CA.GG.TT.GA.GA.AA..T..T.TT.GA..T.TA.GG..T.GG.AC.GC.GC.GT.GA.CC..T..T.GA.AA.AA..T.GC.TA.AA.GG.GG.TT.AC.AA.AC.CA..T..T..T.CC.GA....GT.GC.AT.GA.AA.GG.GG.GT.GA....AA.......AC.GA.CA..G.GG.GC.CC..G.GT..T.GA.GG..T.AA.GC.GT.GA..G.GA.GG.GA.CC.AA.AC.AA.GA.AC.AA.GC..G.AT.GA.AT.GG.GC.GT.AA.TT.GG....TT.TT.GT.AA.AC..T.TT.GA.AC.CC.AT.CC.GA....GA.GC.GG.CC.AA.GT.TA.AA.GA.GG.AT.AT.GA.GT....CA.GA.GA.GG.TA.CA.GT.....T..G.GC.GC.AT..T.AA.GT.AA.GC..T.GC.GG.GT.AA.CA.GA.GG.GG.GC.TT.GA.GC.TA.AT.AT.TT....GA.AC.GA....AC.TT.TA.TA..T.AA.GG.......GC....GA....AC.AC.GT.AC..G.AA.GA.TT.....G.AC.GG.GA..T.GA.AT.TA.GG.AA.GT.AC.AT..G.GG.GC.GG....GG.AC.AC.TA.AT....GG.GG.AA.GC.AC.GC....GT.GC.AA.AA.GA.AA.AA.GG.AA.AA.CA.AC..T.AC.GT.GT.AC....CC.GC..T.AA.GA..G.AT.TT.CA.GT.CA.GC.GG.GC..G..T.CA..T.AA.CA..T....GT.AT.AA.GG.CA.GC.GT.AA.AA.AT.GA.GG.AC.GC.AC.AC.GG.TGGTT.GC.GT.GA.AC....GG.GC.GA.AC.AA.GG.TA.AA.GA.CC.GT.AC.AA.GT.AA.ATGCC.AA.TT....GA..T..G.TT.TT.TA....GA.GA.CA.....T.GT..T.GC.GG.AC.GT.GA....GG.GT.GA.AC.AA.GC....TGGAT.GA.AT.GA.TA.GA.TT.GA.CA.GT.....T.AA....GA..G.TA.TA..T.AA.AA.AT..T.AC.GT.GG....AA.AA.TT..G.AT.GT.AT.GA.AA.GA..G.GT.AA.GA..T.ATGTT.CA.GA.ATGGA.GG.AA..T.GT..T.CC.AA.AC.AT.AA..T..G..T.GC.TA.CC.GG.AT.GA.GC.CA.CC.GT.GT.TA.CA.CA.GA....GG.AA.GG.GG.GG.AT.TT.AA.AA.GG.GG.GG..T.AC..T.GT.GA.GC.GT.AT.GG.AC.AC.GG.GG.GG..T.GC.GG.GA.AA.TT.GA....CC.GG.AA.GA.GC.GA..T.CA.GG.GG.GG.AT.TA....AA.AA.AA.GG.GG.TT.AC.GC.GT.GT.AC..T.AT.AA.....G.AT.AA....AA..G....AA.CA.GG.GC.GG.AT.TA.AT.AA....GG.ATGGC....GT.AC.GA.AA....AT.GT.GG....AA.AT.GT.GC.TA.CA.GG.GG.GG.AT.TA..T.AA.GA.GG.AA..T..T.AT.AC.GA.AA....CA.AT.GC.AA.CA....CC.CA.AC.AT.GA.TT.AC.AA.GG.AT.GG..T.TA..T.GC.GG.GG.AT.GT.AC.AT....GG.GG.......AT.CA.CA.AA.AC....TA..T.AA.......GG.CA.GG.GG.GC.AT.TA.AA.GC.GC.GA..T.AC.AT.AA.AA.GG....CA.GT.GA....AA.AC.AC.GT.ATGAC.AC.GG....GG....TA....GG.AC.GG.AT.TA.AA.GC.GG.CA..T.CA.AT.AA.GA....TA.GT..G.AA.AA.TA..T.TT.GA.GG..G.GG.AC.GG.AT.TT.AA.GC.AA....GG.AC.GT.GT.AT.AA..G....AA.AT....GG.AA.CA.AC.CC.CA.GG.CA.GG.GC.GG.GC.GG.AT.TA.AA.GA..G.GG.AA.AT....AT.CA.GA....AT..T.CA.AA.AA.AA.AT.AC.TT.GG.CC.GG.GG.GG..T.TT.AA.AA....GG..T.GT.AC.AT.GA.AA.AC.GT.TT.GC.AA.AA..G.....T.AC....GG.GG.GG.AT.GT.AA.TT.AA.AT.GA.GC.GT..T.AT.GT.AC.AA.....T.TT.AA.AA.AA.TT.GG.CA..T.CA.GG.GG.GC.AT.TGGAA.AC.AA.GG.AA.GT.AT.AT.GA.AA.......TT.TA.CA.AA.GT.AC.GC.AA.AA.AA....AA.GG.GC.GC.AT..T.CA.TA....GC.CC.AA.CA.AT.AC.TA..T.AC.AC.AA.GC.GG....CA.GC.GA.AA.AT.GT.GT....GA.AT....GC..T.GA.CC.CA..T.CA.AT.AA.AA.GA.AT.AC....TT..G.GT.CC.TT.GA.AT..T.AT.AA.GA.AA.AT..T.AC.GT.AC.CA.AT.AC....GT.GC.GG.GG.TA.AA..T.GC.GT.AA..G.GG..G.AA.CA....AT....GC.AA.CA....AC.GG....AA.GT.TA.GT.AT....GG.....T....AT....AA.AC.AC..T.GC.GC.AA.GT..T.CA........T.GA..T.AA.GG.GT....AA.TT.GG.TA.GT....AA.GG.GG.GG.GT.TT....ATGAC.GA.AT.GT.AC.GC.CC.AC.AT.AT....AA....AC.TT.TA..G.AA.GC.GC.GA.AA.GG.GC....GC.TGGTT.GG.CA..T.GT.AC....GA.CA.ATG......AA.TT.GT.AA.GC..G.CA.AA..........CC.TA..T.AA.AA.AA.AT..T.GC.AA.GG.TA.AA.GC.AA.TT.GC.GG....GC.AA.GT.GA..T.CA.CC.AA..T.TGGAA.GG.GA.GC.CC.TA..T.TA.GA.AT.TA.GG.GA.TT.GT....GG.GG.AA.AA.TT.GT....GC....GC..T.GG.TT.AC.AC.GT.GA.GC.CC.AA.GG.AT.AA.AA....GC.AC.GT.AT.AA..T.GA....AT.......TT.CC.AA.CA.GG.AA.TT.GT.AT..G.AT.TGGGC.CA.GA.GG.TT.GC.TT.AA.AA.ATGAA.GT.AC.GC.AT.AA.AC.GT.CA..T.AC.AA.AA.AC..G.TA..T..T.AC.GT.AA..G.GG.TT.AA.AA.AC.AC.GC.ATGGC.AT....GA.GG.GC..T.GT..T.CA..T.GA.AA.....T.CC..T.AA.CC.AA.GA.TA.TA.GG.GT..T....CA..G..T.GA.CC..T..T.GG.....T.CA.TA.AC.GG.CC........T.....T.AT.CC..T.CC.GG....CC.GT.AT.GG..T.GG.AA....AA....TT.AC.GG.AA.AA.AC.AC..T.AC.GT.GC.AT.AC.GG.GA.TA.GA....GA.AT..T..T.GT....GA.TT.AC.GT.TT.GA.TA.AA..G.GA.AA.GA.GA.AA.GT.AT.CC.TT.CA.AT.CA.GT.GG....AA.ATGATG.T.GT.AC.GG.GT.AA.AA.GT.GC.GG.....T.AC.GT.CA..G.GG..T.TA..G.AC.CA.AC.TA....CA....GC.GG.GC.....T.AT.CA.GG.TT....GG.AT.CA.AC.GT..T.GT.GC.GA.AT....AC.ATGGA.CA.AC..T.AT..G.GT.GG.GC.GC.AC.GG.AT.CC.AC.AA.GG.GT.TT.AA.AT.CA.AT.AA.GA.AA.CA..T..G.GT.AC.CA.GT....TA.GA.AA....GC.GA.AA.TA..T..T.AC.GT.CA..G.GG....AA.GG.AC....GC.GC.AT.CA.GT....GA....AA.GT..T..T.AT.AC.GA.CA.AA.AA....AA..G.AA.AT.GA....AA.GG.GA.GC.AT.GC.CA.ATGGA.AT.GG.GC.GC.AA.GC.GC.TA.TT.CA.GT.AA....AC.GC.GA.AC.AT.GA..T.AA.CC.GG.AA.GG.AT.GT.GA.AC.GC..T....GG.AA.TT.....T..G.GC.GC.AT.ATGAA.GC.AA.GC..T.AA.GG.CC.GC.AC.AT.TA..T.CC.GC.GG.AA.TA.AC..T.GA..T....GG.AT.GA.CA.GG.GA.TT.GC.GG.GA..T.GA.AT..T....AA..T.AC.AT.AT.GG.GA.GC.GC.AA.AC.AC....AT.AT.GC.AA.AA.AC.GA..G.AT.TT.GA.AT.TT.CC.GG.GT.AA..T.AC..T.AA.GG.GT.AC..T.AA.AC.GG.GT.GG.GT.AA.CA.....T.GA.GG.GA.GG.GG.GC.AT..G.GT.AA.CA.GG.AA..T.....T.AT.AA.AA.ATGAT.AC.GG.AC.CA.GC.AA.CA.GG.GG.GC.AT.TGGGC.AA.GG.GC.AC.GT..G.AT.GA.AA....AC.AT.TA....AC..G.GC.AC....AC.GG.GG.GC.AT.GC....TA.AA....GA.GT.AT.ATG...AA....AC.GT....GG.AA.GG.GT....AC.GT.AC.CA.GG.GG.GG.AT.TA.AA.GG.GC.TT....AC..T.AA.GT..........AC.AT.GT.AA.AA.GA.GG.GC.AA.GG.GC.GG.AT..T.AA.TT.GG....GT.AA..T.AT.AA.AC..T.GT.GC..T.AA.CA.GG.GA..G.GA.GT....GG.CA.TT.GA....GA.GG.TA.AA..T.AT.GG.AA.GC.AA.GG.GC.AC.GG.TT.......GG..T.GG.GA.CA.TT.GG.AC....GG.GC.ATGAT.GA.CC.GT.AT........T.GG.GA.CA.GG.GG.AC.AC..T.AC.CA.GC..T.AA.GT.GG....CC.GC.AT.GA.GC.GG....GG..T....AA.TA..G.GG.AT..T..T.AA..T.AA.....T.CA.AT.GG.AA.AA.AA.AT.TT.AT.CA.GT..G.GA.AC....GG..T.CC.GC.GT.CC.TT.AA.AT..T.GT.GA.GC.AC.CC.GA.AA.CC.GA....GG.GA.AA.AA.AT.ATGAC.GT.AC.AA.AT.AA.GA.AC.GT..T.AC.GT.AA..G.GG.GT.AA....AC....GC.AA.AC.CA.GT.AA.CC.GT....TT.CA..T..T.GG..T.CA.GA.CC.AA..T....GG.AT.AA.GT..T.GT..T.AT.GA.CA..T.GG.AA.CC..G.TA.CA.AA..T.GG.ATGTA.CA.AC.GA.GC..T.GA.AT.GG.GC.TA.AA.GG....GT.TT.GC.AC.CC.CC.TT....AT....AA.CC.AT.AT.AA.AT.CC.GA.TA.AA.GG.AA.GA.AT.CC.AA.GG.AC.GC.AA.GG.GC.TA.AT.GA.GT.AT.TT.AA.AT.AC..G.AC.GA.CC.AA.TA.GC.CC.AA.AA.AT.GG.TA.TA.AT.AA.CC.GC.GG.AA.CA....GC.GG.CC....GA.TT.TT.AA..G.....T.GG.CA.GT.TA.AA....TA.CC.AC.GG.AA..T....TT.AA.....G.AA.GA..G.AA..T.GT.TT.CA.AT.AC....GA.GC..T.CC.AA..T.GA.AA....TT..G.AT.GT..T.TA.GA..T.AA.AC..T.AT.CC..T.AC.GA.GT.......AA.GC.AC.AT.CA.GA.GA.GA.GA.GA.CC.AT.AT.GT....AT.GT....AA..T.AC.CA....AA.GA.GG.AC.TT.GG.AT.AA.GT.AC..T.GA.CA.GC.GC.GG.AA.AA.AT..G.GT..T.GC....AC.GA.AC.GT.GG..T.AC.GC.AC.GG.GG.GC.GA.TT.AC.GT..T.GA.AA.TA.AA.AT.GT.TT.GC.CC.GG.AA.AC....CA.TT.AT.CC.AT.AA.GT.TT.AA.GA....AC.GT.AA.AC.GA.AA.CA.TT.AA.GT.AC.AT....GG..T.CA....AT..T..T.CC.AC.AA.GA.AT....AT.GG.GG.GG.GG.GA.AA.TG..T.GT.AC.AT.AC.AA.GA.GA.AA.GC.....T.GT.GT....AA.GT.AC.AA..T.AA.GG.GA.AA.GG.AA.GA.GC.CA..T.GT.TT.CA.GT.AC..T.AC.AA....GT.GC.GA.AC.GT.AC.GT.GA.TT.AC.AC.GT.AA.CA.AC.GC.GT....GG.AA.GA.TA.GT.CC.CA....GG.AC..T.GT.TT.AC.GC.AA.GG.GA.CA.AC.AT.CA.AT.GT.GT.CC.GT.AT.AA.GA..T..T..T.AA.GG.GA.AA....TT.TT.GT.CA..T....AA..T....AC.GG.GT..T.......AA.AA.TT.....T.CA.GA.GG.GT.GG.AC.AT..T.GA.GA.GA.GC.CC..G..T.AA.GT.GT.GA.GT.AC.AA....AA.GA.GG.AC.TT.AC.TT.AA.GT.AA..T.AC.AA.TA.AC.GG.GC.GT.AC.GT.GA.GT....AC.AT.GA.GG.AC.GC.AT.GC.GG.AA.GG.AA.GA.TA.TA..T.AA.AA.AA.GG.GT.GC.GT.GA....AT.AC..T....TT.AA.GT.AA.CC.AA.AT.GA..T.AC.GG.AA..G.AA.CA....TT.AC.GT.CA.GT.AA.GA.GA.AC.GA.TT.AA.TA.GA.AA.TT.TT.GG.AT.GT..T....AA..T.AT.GC.AA.GG.AA.AA.AA.GA.TT.AT..T.AA.GA.GG.AT.GG.AC.AT.CA.AA.GA.GA....CC....GT.AC.AT.AA..G.GT....GG....AA..G....GG.TT..T.CA.TT.AC.GC.AT..T.GC.GG.GA.AT.CC....GG.TT....AC.GT.GT.AA.GC..T.GA.CA.GA.GC.GT.AA.GG.GT.GA.TA.AT.TT.GA.CC..T.AC..T.GA.TT.GA.GG.TA.AA.GG.AA.AC.AT.AC.TT.CA.GT.AA..T.AT.AA.GA..G.AT.AA.AA.TGGAC.CC.AC.GC.AA.TA.AA.AA....TT.AA..T.AA....AA.AA..T.AC.AT.GG....GA.AC.AT..G.GG.AT.GC.GG.AA.GG.AC.AT.TA.GA.GA.GA.TT.TA.GT....GC.GG.TA.GA.GT.GA.AT.AC.TA.AT.CA.AC.GG.GG.AA....AA....GA.AA..T.CA.GA..T.CA.CA..T.GT.GC....AA.TT.AC.AA..T.GG.AT..T.GC.CA.CA.GA.CC....GA....AA.AA.GG.AC.GA..T.GC.GT.CA.ATGAA.GC..T..G.AA.AA.AT..T.AA.GT..T..G.CA.AT.GA.CA.TA.AA.GC.GA..T.GC.GA....CA.CA.GT.CC.AA....GT.GT.AC.GT..G.AT.GA.GA.GT.AC.CA....AA.AA.GG.AC.TT.AA.TT.AC.GT.AC..T.GA.AA.GC.GC.AA.AA.GA.AT....GT..T.GC....AT.AA.TT.......GC.AC.GC.GC.GA.TT.GC.AC.CA.AA....CA..G.AT.AT.AT.CC.GC.GG.AA.AC....GC.AC.TT.AC.AT.GA.GT..T.GA.GA.AC.AA.GT.AA.GA.AC.AA.GT.TT.AA.GT.GC..T....GA.GC.TT.TA.GG.GC.GA.AA.AA.....G.CA.GT.CA.AT....GT.GA.GC.GG.GT.GG.AC.GG.AC.AT.AA.GA.AA.GA.GT.GC.GG.GT....AT....GA.AC.AC.GT.AT.AA.GG....AC.GC....GT.AT.AT.....T.AC.AA.CC.GT.....G.GA.GT.AC.GT.AA..T....GT....CC.GA....GC....GC.GC.GA.TA.AA.AC.GT.CC....AA.GT.AT.AT.CC.GC.GG.CA.AC....GT....TT.AA.AT.AA.AC.AA..T.GA.AA.GA.GT.AA.AA.GA....GG.GA.AA....TT.CA.AT....AT.AC.AA.GC.TA.AT.GA.GA.GG.CA.GG..G.AA.GA.GT.AA.AT.AA.AA.GG....GC.GT.AT.AC.AT....GA.TT.AA..T.AA..T....AT.......GT.AC.GC.AA.AA.GG.GG.GT.TT.AC.TT.GT.GT.....T.AA.CA.CA.GT.CC.AA....AA.AC.AT....GT.AA.GT.AA.AC.AA....GG.GG.GA.TT.GG.AC.GC.AC.GT.GG.AA.GA.TA.AC.GC..T....AA.AA..T.GT.TT.AA.GC.GG.AA.CA....AA.AC.GT.AC.GT.AT.GT.AA.GA.GA.AA.GT.GT.AA.GG.GA.AA.AC.TT.AC.GT.AA..T....AA.GC.AA..T.AC.GG.GA..T.AA.CC.AA..T..T....GA.AA.TA.......GC.AC.CC.AC.AT..T.AA.GA.GA.AA.TA.AT.GT....AT....AA.AA..G.AT.AA.AA.GA.CA.AC.ATGAC.TT.AC.AT.....T.GA.CA.GC.GC.GA.AA.GA.GT.AT..T.TA.GC....AC.CA.GC.AC.GG.AC.GG.TA.GG.TT.GC.GC..T.AA.AA.GA.TT.GA.CC.AT.GT.AA.AA.AA.AT.AC.AT.GT.AA.GG.GA.AC....GT.AA..T.CC.GT.CC..T.GG..T.AA.GA.TT.GC.AA.GG.GA.AA.AC.TT.CA..T.GT.AT.AC.GG.GC.AA.TA.GG.AC.AA.AA.GG.GC.GC..T.GG.GC.CA.TA.GT.GA..T....AA.GC.GT.GG.AT.GG.AC.AT.AA.GA.TT.GT.CA.AC.AT..G.AT....GA.GC.AC.GC.CA.AA.GG.GG.GT..T.AC.TT.AC.GT.....T.GA..G.GC.....T.GT.CC.AT..T.GT.AA.GC.GA..T.GT....GG.TA.GC..........AA.TA.GG.AC.AT.AA.GG.GT.AA.AA.GC.GG.ATGAT.AA.TT.GC.GT.GA.CA.TA.AC.GG.AA.AC.GC..T.AC.AA.AC.AT.AC.GT.AA.AC.TT.GA.GA.AA.AT.GT.AA.AA.GA.AA.AC.TT.GG.AT.AA..T....AA.GC.CA.TA.GA.AA....GT.AA..T.GT....GC..T.AA.AT.AA.GA....GT.GG.AT.GG.AC.AT.AT.AA.AA.GA....GT.GC.AT.AC..T.AA.AA.GG.GT.CA.AA.AA.GA.GG.AC.GT....TT.AA.AT.AA..T.GA.AA.GC.GC....AA....AT.AC.GT.AT.GC.AA.AC....GG....GA..T.AT.AA.GC..T..T.GA..G.AA.AT.GT.AT.CC.GC.GG.GC.AT....AA.GT.GT.AC.AT.....T....AA.AC.AA.....T.CC.AA.GA.GT.AC.TT.AC..T.AA.AT.......CC.CA.TT.AA.GG.CA.GA.GC.GG.GA.TT..T....AT.AC.AA.GA....GC.AC.GG.AC.TT..T.AA.GG.AC.GT....AC.GG.GG.GA.AA.GG.AA....GT.AC.AC....AT.AC.GG.AT.AC.GA....GA.TT.CA.GT....CA.GC.GG.AT.AC.AC.TT.AC.CC.TT.GG.AA.GC.GG.GC.......GC.GA.AA.TA.TA..T.GC.TT....GG.GC.GC....GG.AC.AC.....T.GC..T.AA.TA.GG.AC..T.AC.AT.AA.CA.GA.AA..T.AC.GT.AC.TA.AC.CC.GA.GA.AA..T.AC.....G.AC.AC.GG.AA.AT..T.CA.CC.AA.TA.GA.GG.AA.GA.AA..T.AT..G.TT.AT.GG.AA.AA.GC.TT.GT..T.....T.TT.AT.AA.AA....GA.GG....GG..T.CC.GC.TT.GT..G.GA..G.....T.AC..T.GC.GT.AC.AA.AC..T.CC.GT.TT..T.GG.....G..T.CA.GA.GA..T.GT.AA.AA....GA.CA.GA..T.GC.AC.GC.AT.CA.AA.GC.AT.AA.AT.TA....GG.GG.GA.TA.AT.GT.GA.AT....AA.TA.TT.TA.GA.CC.GA.GG.GA.AA.AT.AC.AT.AA.GG.AT..G..T.TT.GG.GC.CC.AA.GC.GA.GC.GG.AC.GC.ATGAA.GT.AA..T.GG.AA.TGG...CA.GA.GG.AC.TA.TT.AA..G.GC..G.TA..T.AA.GT.AA.GT.AT.GA.AA....AA..T..G.TT.TT.AA.GA..T.AA.AA.AC.GA.AA.GA.AC.AA.GG.AC.CC.GG.AT.TA....GC.AC.GC..T.AC.AA....TT..T.AT.AC.GT....GA.GG.CA.AC.GT.AC.GA.GG.GT..T.AA.AA.AA.GG.CA.CA.GA..T..T.AT....AC..G....GT.GG.AA.CC.AT.CC.AA.AC.CA..T.GA.CA....TT.GC.GT.AC.GT.AC..G.CA..T.AC.CA.GA.AA.AA.GC.AC.AA.GC.AA.AA.TT.AC.GT.GA.AA....AA.TGG....T.GA.GC.GC.GC..G.AC.....G.AA.GG.......GA.TA.CA.GA.GT.GG.GG.GT....GT.GA..T.AA..G.GG.AC.TT..T.GT.GG.CA.GA..T.GC..T.GA.GG.GC.GG.GT.CC....AA.AA.ATGAT.CC.GG..T.GT.TA.GA....GG.AC.GT.AA.CC..G.CC.GT.AT.CA.GC.GT.GT.AA..G.AA.GC.GA.CA.CC.GT.CC.GT....AT.AC.GC.AC..T....TGGAT.GA.CA.CA.AA.AA.CA.GA.AA.GT.GT.GT....AC.AA.AC.TA.AC.AT.GC.GC.AA....GT....AA.GA.AA.GA....AA.TA.AC.AT.GC.AT.CA.CC.GA.AC....CC..T.GT.AC.GG.GT.TA.GG.TGGAA..T.AA..T.AA.TT.AA.TT.GC.GA.GC.AA.AA.AA.AA..T.AA..T.CA.CA....GG.GC.GT.CC.GT.GT.GT.AA.AA.AC.GT.AA.CA.AC.GA.TA.GC....AA.GG..G.GT.AA.CA.CC.AC.GA.GC.AT.TT.GG.AA.GG.TGG....T.GC.GG.AT.CC.AC..T.TT.AT.GA.AA.CA.GG.AA.AC.CA.AC.GC.GA.GA..G.AT.AT..T....TT.CC.GG....AC.CC.CA.GT.TT..G.CC.AC.AA.GC.AA.CC.AA..G.TT.GA.TA.GG.AA.CA.CA.TT.AC....GT.GA.GT.TA.....G.TT.GC.GC.AT.GA.CA.GT.AA.GG.TT.GG.GC..........CA.TGG...AA..T.GG.TT.GG....GA....GA.CC.GG.AA.TT.GG.AC..T.AA.AC.AT.GA.GG.GG.AA.AA.AT.GT.TA.AC.GA.GT....GG.AC..G.TA.AC.TT.AC..G.TA.AA.AT.AC.......GT.AA.GA.TA.AA.GT.TA..G.AT.AC....GT.AA....TA.GT.CA.GC.AA..........GT.AA.GG..T.GT.TT.GG..G.GA.GC.AA.GG..T..T.AA.AC.AT.AA.GC....GA.AA.AC.AC.TT.AA..T.GA.TA.GG.GG.AC.GA.CA.GT.CA.AC.AT.AC.GC.TT.GA.GG..G.GT.AT.AA.TT.AA.TA.GC.GA.CA.AA..T.AA.CA.AT.AA.AA.GT.AA.GG.TT.AA....GG.GC.CC.GC.....T........T.GC.GA.GC.CA..T.AC.GT.CC.GA..G.AT.CA.CA.TT.AA.TA..T....GA.GG..T..T.AC.AA.AA.AA.TGGTA....GG....GT.AC.GA.AT.AA..T.AA..G.AC.AC..G.TA.AA.TA.AA.AC....AC.AC.CA.GA.AA.GA.AA.CA.GA.AT.GT....GG.AC..T.AA.AA.AT.AT.TA.GG.AA....AA.GA.GG.AA.AT.AT.TA.GA.AT.AA.GC.GT.GC.CA.GC.GG..T.AA.AA.AA.GC.AC.AC..T.AT.AC.AA.GA.GA..T.AT.GC....GT.AC.GT.AA.AA....GA..T.AA.AA.TA.GT.GC.GA.AC.AA....GG.AC.CC.AA.AC.GG.AA....AT.AC.AA.TA.AA.TT.GA.AC.AA.GG.CA..T.GT..G..G.AA.AA..T.TT.AA.AA.AC.....T..G.AA.AA.AA.TGG...TA.GA.CA.GT.GG....GT.AA.AC.TA.AC.GA.CC.TT.GG.CA.ATGAC.TA.TA....TA.GA.TA.AA.AT.CC.GC.TT.TA.AA.AA....GA....GA.AA....GA.AA.AA.GA.CC....CC..G.TA.GA..G.GA.GA.TA..G.GG.AA.GT.GT.AA.AT.AA.AA.AA.GA.GG.AA.AC.AT.TT.AA.TA.AA.AC.GA.GA.AA.AA.AA.AA.GC.GT.GG.AC..T..T..G.AC.GT.AA.TA.GC";//.GG.GT.AA.AC.GT.TA.AC..G.GA....AA.GG.ATG.T.AC....AC....AT.CA.GG.GA.AA.TA.AC.TA.AC.AA.AA.TA.AA.TA.TA.GA....AC.GA.GG.GT..T.AA.GG..T..T.AA.AA.TA.AC.GA.AC..T.GG..T.AA.AC.TA.TA.AT.TA.AC.AC.AT.AA.GG.GT.ATGGA..G.CA.GT.GA.AA.AT.AA.AC.GT.GA.AC.GT.AC.GG.CA.AC.GA....AC.AA.TA....TA.GA.AC.TA.GG.AA..T.GC.AA.GT.AA.AC.TA.TA.GG.GA.GA.GC.TA..T..T....AC.AC.TT.AA.AC.TA.GA.AA.AC.GG.AA..T..G.GA.AC.AT.AT.AA.GA..G.AC....AA.AT..T....CA.CA.CA.TT....TA.GA.GC.TA.GG.GC..T.AC.GC.AT.GT.AA.GC.GA.GG.AT.GC.AC.AT.AA.GC.TA.GA.AC..G.GG..T..T.AA.CA.AC.AT.GA.GC.AA.GG.AA.AC.TA.AA....CA.TA..T.AA.AC.....T.GC.AT.CA.CA.GT.AC.AA.TA.AA.TA.TA.GC....GG.GC..T.AA....GT.GC.GT.AA.AA.AA.AA....AA.GG.AC.GT.CA.AA.TT.TA.TA.G";//A.CC.AA.AA.CA.CA.CA.TGGAT.TA.AC.AA.AA.AA.CA.GC.CA.GC.CA.GT.AC.CA.AT.AC.AC.GA.AA.TT.GG..G.AC.GA.AT.AA.AA.GA..T..T.AC.GG..T.AA.AC.AA.TA.CA.TA.GA.AA.CA.CA....AA..T.GT.AC.CA.ATGAA.AC.CA.GT.AA..T.AA.AA.GA.GG....AA.GG....AT.GA.ATG.G.TA..G.TA.GA.AA..G.GG.AA..T.AT.TA.AA.AA.CA.CA.GG.GC.AA.AC.AC.GC....AT.AA.AC.TA.GA.AA..T.GG.AA..T.AT.AA....GG..T.GC.GA.CC.CA.GG.GC.AA.TA.AA.TA.GT.AC.AA....GC.GG.CA.GT.AC.AA.CA.AT.GG.CA..G.TGGAC.TA.GA.GC.AC.AC.AA....TGG...AA.AA.AT.TA.AC.AC....TA.AA.TA.GA.CA..G.GG";//..G..T.GT.AA.GT.GT.GA.CC.AT.GC.GC.AA....GA.AC.GA.AT.GA.TA.CC.GC.AC.TT.GT.TT.GA.TA.GA.AC.GC.AA.CA.AA..T..G.GT.AC....AC.....G.AC.GG.AT.GT.AC.AT.CA.TGGTT.GA.GC..T.GG..G.GT.GT.AA.GA.CA.AA.GC.TA.GG.GG.AT.AC.AC.TA.AC.TA.GA.GA.GT.GG.AA..T.AA.CA.AA..G.TT.AC.CC.GC.CA.AA.GA.CA.GA.GC.AC.GT..G.AA.AC.AA.TA.CA.TA.GA.GG..T.GG..G..T..G.GC.GT.AA.TA.CA.AA..T..T.GG.AT.GC.CA.TA.TA....GC.AT.GA.TA.TT.AA.CC.GG.CA.AA.AA.AC.AC.GT.GA.AA.GT.CA.CA.TA.TGGGT....GC....TT.AA.GC..T.AC.CC.AC.AA.TA.GC....AC.GC.AA.CC.AC.AA.GC.AC..T.AC.GG.CA.AA.GC.AA.GG....CA.GT..T.CA....CC.AA.GG.AA.CC.GA.GA.CC.GG.GA..T.AT.AC.AC.AT.AA.AC.AT.TA.AA....GA....AC.TA.CA.GT.AA.AT.AC.AC.CC....GG.GA.GT.AA....CA.AA.AC..G.AC.AC.AC.GC..T.TA.AA.GC..T.GG.TGGGT..T....AC.AC.GA.AC.GA.GG.GA.CA..T.GT.AC.AA.AA.TA.AA.AC..T.GG.CA.GT.GT.CA....AC.AA.GG..T.GG.CA.GT....AC.TA.AC.TA.AC.....T.CA.GG..T.GT.AA.AC.GT.AC.AC.GT.GC.GG.AC....AA.TT.AA.TA.GA....GC.GG.AA..G.AC.TA.CA.AA.AA.GT.GA....GA.AT.AA.AC.TGGAC.TA.GA.AT..T.AA..G.GT..T..T.AA.AT.GC.AC.AT.GA.AC.TT.AA....TA.GG....CC.AA.GG.CC.GG.CA.GT.GT..G....TGGAA.TA.CA.GG.AA.GA.AT.GT.TA.GT.AA.TT.AA.GG.AA.AA.GT.TA.TA.GA.AT.GA..T.AT.CA..G.CA.AA..G.AA.AT.TGGGA.GG..T..T.TA.GC.AC.TA..T.TT.AA....GC.GG.AA.AT.AT....GC.GA.AC.GC.GA.AA.GT....GA.AT.CA..T.AC.TA.AA.GA..T.GG..G.GT.AT.AA..G.AA.AC.GT..T.TT....AA.AA.GT.AC.AA.GA.AT.GG.TT.GA.TA.AA.GG....AC.AT..G.AA.AC.TA.CA.TA.AA.GA....GG.AT..G.AA.AA.TA.AC.AT.AA..T.GA.GC.AC.GA.AA.AT.GG.TT.AC.AA....TA.GT.AT.GA.AA.GG.AA..G.AT.AT....ATGAC.CA.......TA.ATGAA.CC....AA.AA.AC.CA.AA....AA.AA.GC.GT....TT.AA.TA.AC..G.GA.GG.CA.AT.CA.AA.GT.AA.GG.AC.CA.GA.AT.GT.AA.AA.AT....GG.CC.AA.GC....AT..T.AC.AA.TA.GC.TA.AA.AA.GA.GG.AA..T....GC.GT.AT.AA....TA.GC.TT.AT.AA.GC.GG.AC.AC.GC.CC.TT.AC.GT....AA.AT....GA.TA.TT..T.CA.TA..G.AA.GG.CA.GC.GA....AC....CA.AC.AT.CA.TT.GA.AA.CA..T....GG.GT.AA.AC....TT.TT.GG.GG.....T..T.AC.AC....TT..G.TT.GA.GA.CA.AA.CA..T.AA.AC.GT....CC.AA.CA.....T.......GC.AC.AA.AC.GC.TT.AA.TA.AA.AA.GG.AA.AT.GT.TA....GC.AA.AC.AA.AT.GG..T..G.GA..G.GT.CA....GA.AA..T.TA.GA.TA.AA.TA.GA.AA.AA.GG..G.AT.AT.AA..G.AT.TT..G....CA.TA.GA.AT.GT..T.GA.GC.AA.GT....AA..........AA.TGGGT.AA.GC.AA.TT.AC..T.AA.GA.AC.TT.GA.AA....AT....AA.GT.AA.GA.GC.CA.GC.GT.AA.GG....CA.AT.GT.TGGTA.GA.AA..T.AC.GT.AA.GA.AC.GC.AC.TT.AC.AT.AA.GA..T.GA....GG.AA.TA.AA.AT.TA.GT.AC.TGGGT.CC.AT....AA.GG.GC.AC.GC.GG.AA.TA.TT.AC.GT.AA.GG.AA.AA.AT.AA.AA.TT....AT.GA.TT.AC.AA.CC.GT.TA.CA.AA..T..T.CC.AA.AA.CC.GA.TA.AA.AA.CA.CA.TT..T.AC.AA.GA.GG.GT.GA.TGGAT.TG..T.AA.AC.GT..G.CC.GA.......GG.AA.AT.AC.GT.CA.GT.GG.GG.AC....TT.TA.GG.GG.GC.TA..T.GC..T.GA.GC.AT..G..T..G..G.AT.GA.AT.GT..G.GA.GT.TT.GC.GC.GG.AA.AC.AA.CA.AC.AT.AA.TT.AA.......AA.CC.CA.AA.GT.AA..T.AC.GG.AA..G.TGGGA.GA.TT.TA.AA....GG....AA.CA.TGGTA..G..G....AC.AC.GG.GG....AC....AC.AA....GT.AA.CC.GT.GC.AC.TA.AC.TT..G.TA.AT.CC.GG.GG.CA.TA.CA.AT.TT.GC.GA.TGGCA.GG.GG.GT.AA.AA.GG.GC.AT.TA.GT.AT.AA.AA.GG.GG.CA.AT..T.GC....GT.TA.GC.GA.CA....AC.GG.AA.GG.AT.TT.AA.AA.CA.TT.AT.GC.AA.AT.AC....AC....GG.TGG.T.AC.AT.AA.AT.CA.GT.CC.GG.GG.AA.CA.GC.AT.GC....AA.AT.AA.AT.GT..G..T.GA.AC.CC.CC.GT.AC.AA.TA.GA.TGGGA.TA..G.GG..G..T.AC....GT.ATGTA.AA.AA....GA.GC.GG.CC.GT.AT.AA.CA.GT.AA.TA.TT.TA.AA.GC.TT.GA.CA..G.GT.GG.CA.CA.GT.AA.AC..T.AT.GA.GA.AA..T.GT.CA.AC....TT.AA....TA.AC.TA.GA.GG.AA.AA..T.GT..T.CA..T.GG.GA.AA.AC.AA.GG.TT..T.AT....CC.TT.GG.CA.AC..T.GC.GT.GA.....T.AT.AA.CC.AA....AA.GC.GA.GG.GA.AC.....T.AA.AC.TA.TGGGC.GT.TA.GA.TA.CA.GG....GC.AT.TT.TGGTA.CC..T.TA.TT.AA..T.GA.GG.CC.AT.CC.GA.AA..T.GT.TA....CC.AC.CC.GA.AA.TA.TT.TA....AC.CA.GG.GT.GT....AT.GT.AA.CA.GG.GA....CA.TT.AT..T.CC.GC.AT.TT..G.GG.AA.CA.TA.GA.GA.AA.AC.GG.ATGTA.TT.GT.GG.AA..G.GC.TA.GA.CC.CA.TT.AA..G.TA..T.AC.AA.CA.GT.GC.GG..T.GA.GC.GG.AA..T.TA....TT.GC.GT.AA.GC.CC.GT.AC.CC.AA.AT.CC.GG.GT.CC.AA....TA.AC.AC.AA.GA.CC.GT.AA.AA.CC.GT..T.......TT....TA.AA.AC.GT.GG.TGG.T.GT.GG.GA.GA.GT.TA.TT.AA.GG....TA.AC.AA.AT.GC..T..G.GC.GC..T.AC.CC..T.AA.GT.GT.GG..T.GT.GA.CC.TT.GG.AT....GA.TT..T.GT.GG.AC.GC.TA..T.GC..T.GG.AA.AA.AA..T.GC.GC.AT....TA.GC.GG.AT.GC.AT.CC.TT.GG.TT.AA.AA..T.GG.AA.AT.GC.GT.AT.GC.GC..G.GC.GC.GA.CA.AT.CC.GT....CA..T.GC.AA....AT.GA.AA.GT.CC.AC.AC.CA.GT.GT.AA.GC.AC....CA.CC.AA..T.GT.GA.CA.AT..T.AA.CC..T.AA.TG.TT.GT.GC.AA.AC.CC.GT.ATG.G.....T.AA....GA.AC..G..T..T.GT.AA.GC.CC.GT.GC.AC....AA.AA.AA.AA.TGGAC.AC.AT.TGG.T.GT.GG.GG.AT.GC.GG..T.AT.AC..T.AA.CA..T.AC.AA..G.AA.AA.AA.CA.TT.AA..G.AA.GC.AC.CC..G.GC.TT.AA.CC.GA.GC.TGG.T.AA.CA.AA.GA.TGGAC..G..T.AT.CC.GT.CC.AA.AA.AC.CA.CA.AA.GC.AA.TGGGA..T....GC.GA.GA.AT.GA.GG.AT.TG.CA.AA..T..T.ATGGG....AA..G..G.......TT.GC.GC.AC..G.CC..G.CA..T.CC.GT..G.AA.AC.GC..T.CC.GC..G.AA.AA.AC.AA.CC.GT..G.AA.GG.AA..G.GC.AT.GT.CC.AA.CA.AA.CA.CC.GT.AT.GC.GG.GA....AA.CA.GT.AA.AA.AA.....T..T.AA....GT..T.CA.AA.";
//		new ORF ().getLevenshtein (a, b);
//	}
	
//	public static void main (String [] args) {
//		OrfRecognizer reco = new OrfRecognizer ();
////		reco.load ("workspace/PALESB58.fasta");
//		reco.loadFromCodonProject (new File ("workspace/PALESB58/CODON.codon"));
////		reco.computeOrfFreeAndOrfSuperpositionSpans (10);
////		for (int i = 0; i < 6; i++) {
////			List <ORF> orfs = reco.orfsPorFrame [i];
////			for (ORF orf: orfs) {
////				if (orf.id.equals ("F1-001107-3168-4275")) {
////					int ind = reco.getCodonIndex (reco.freeSpans, orf.start - 1);
////					int startSpan = reco.freeSpans.get (ind).index + 1;
////					int stopSpan = reco.freeSpans.get (ind + 1).index - 1;
////					if (stopSpan == orf.start) {
////						orf.selectBetterEntryFulfill (orf.lenght + stopSpan - startSpan + 10);
////						System.out.println (orf.lenght + stopSpan - startSpan + 10);
////						System.out.println (orf.lenght);
////					}
////				}
////			}
////		}
//		
////		
////		for (int i = 0; i < 6; i++) reco.orfsPorFrame [i].clear ();
////		ORF orf = new ORF (reco, 4715705, 4718489, 0);
////		orf.id = "F0-002784-4715705-4718489";
////		reco.orfsPorFrame [0].add (orf);
////		reco.loadBaseSequence ();
////		orf.loadProduct ();
////		
//////		System.out.println (orf.start);
////		System.out.println (orf.lenght);
////		int [] bounds = orf.getQueryMatchBounds ();
////		System.out.println ("startquery " + bounds [0]);
////		System.out.println ("  endquery " + bounds [1]);
////		System.out.println ("   lgquery " + (bounds [1] - bounds [0]));
////		System.out.println ("startmatch " + bounds [2]);
////		System.out.println ("  endmatch " + bounds [3]);
////		System.out.println ("   lgmatch " + (bounds [3] - bounds [2]));
////		orf.autoRebound ();
//////		System.out.println (orf.start);
////		System.out.println (orf.lenght);
////		bounds = orf.getQueryMatchBounds ();
////		int rebound = (orf.start - orf.startOriginal) / 3; 
////		int reduction = rebound - orf.product.getStartQuerySeq ();
////		
////		System.out.println ("rebound " + rebound);
////		System.out.println ("reduction " + reduction);
////		
////		System.out.println ("startquery " + bounds [0]);
////		System.out.println ("  endquery " + bounds [1]);
////		System.out.println ("   lgquery " + (bounds [1] - bounds [0]));
////		System.out.println ("startmatch " + bounds [2]);
////		System.out.println ("  endmatch " + bounds [3]);
////		System.out.println ("   lgmatch " + (bounds [3] - bounds [2]));
//	}

}