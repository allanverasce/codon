package codon.data;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class Overlap {
	
	ArrayList <ORF> orfs = new ArrayList <> ();
	public int begin = 0;
	public int end = 0;
	public static int tolerance = 10;
	
	public Overlap (int begin, int end) {
		this.begin = begin;
		this.end = end;
	}
	
	public int getLenght () {
		return end - begin;
	}
	
	public void addOrf (ORF orf) {
		orfs.add (orf);
	}
	
	public void autoCuration () {
		if (orfs.size () < 2) {
			String orfS = orfs.size () + " orfs\n";
			for (ORF orf: orfs) {
				orfS += orf.id + "\n";
			}
			throw new RuntimeException ("BUG " + orfS);
		}
		if (orfs.size () != 2) {
			return;
		}
		ORF orf0 = orfs.get (0);
		ORF orf1 = orfs.get (1);
		if (orf0.start > orf1.start ) {
			return;
		}
		try {
			File f0 = new File (CodonConfiguration.workspace + orf0.reco.dataDestDirS + "/orfs/" + orf0.id + (orf0.hasSwissUpdate () ? ".swiss" : ".unip"));
			File f1 = new File (CodonConfiguration.workspace + orf1.reco.dataDestDirS + "/orfs/" + orf1.id + (orf1.hasSwissUpdate () ? ".swiss" : ".unip"));
			
			// orf0 forward orf1 forward: try to rebound orf1
			if (orf0.isForward () && orf1.isForward () && f1.exists ()) {
				List <ProductProtein> products = ProductProtein.loadAll (f1);
				int lg = orf1.length;
				int lgRef = lg - (end - begin) + tolerance;
				float acc = orf1.getAccuracy () - CodonConfiguration.tolerance;
				ORF selec = null;
				for (ProductProtein p: products) {
					ORF o = new ORF (orf1.reco, orf1.startOriginal, orf1.stopOriginal, orf1.frame);
					o.product = p;
					o.setLoaded (true);
					o.autoRebound ();
					float nacc = o.getAccuracy ();
					if (nacc > acc && o.length <= lgRef) {
						if (selec == null) selec = o;
						else if (selec.isFragment () && o.isFragment ()) selec = o;
						else if (!selec.isGene () && o.isGene ()) selec = o;
						else if (selec.isUncharacterized () && !o.isUncharacterized ()) selec = o;
						else if (selec.length < o.length) selec = o;
					}
				}
				if (selec != null) {
					orf1.product = selec.product;
					orf1.autoRebound ();
				}
			}
			// orf0 reverse orf1 reverse: try to rebound orf0
			if (orf0.isReverse () && orf1.isReverse () && f0.exists ()) {
				List <ProductProtein> products = ProductProtein.loadAll (f0);
				int lg = orf0.length;
				int lgRef = lg - (end - begin) + tolerance;
				float acc = orf0.getAccuracy () - CodonConfiguration.tolerance;
				ORF selec = null;
				for (ProductProtein p: products) {
					ORF o = new ORF (orf0.reco, orf0.startOriginal, orf0.stopOriginal, orf0.frame);
					o.product = p;
					o.setLoaded (true);
					o.autoRebound ();
					float nacc = o.getAccuracy ();
					
					if (nacc > acc && o.length <= lgRef) {
						if (selec == null) selec = o;
						else if (selec.isFragment () && !o.isFragment ()) selec = o;
						else if (!selec.isGene () && o.isGene ()) selec = o;
						else if (selec.isUncharacterized () && !o.isUncharacterized ()) selec = o;
						else if (selec.length < o.length) selec = o;
					}
				}
				if (selec != null) {
					orf0.product = selec.product;
					orf0.autoRebound ();
				}
			}
			// orf0 reverse orf1 forward: try to rebound orf0 &&/|| orf1
			if (orf0.isReverse () && orf1.isForward () && f0.exists () && f1.exists ()) {
				List <ProductProtein> products0 = ProductProtein.loadAll (f0);
				List <ProductProtein> products1 = ProductProtein.loadAll (f1);
				int lg = orf0.length + orf1.length;
				int lgRef = lg - (end - begin) + 2 * tolerance;
				float acc = orf0.getAccuracy () * orf1.getAccuracy () - CodonConfiguration.tolerance * CodonConfiguration.tolerance;
				ORF selec0 = null;
				ORF selec1 = null;
				for (ProductProtein p0: products0) {
					for (ProductProtein p1: products1) {
						ORF o0 = new ORF (orf0.reco, orf0.startOriginal, orf0.stopOriginal, orf0.frame);
						o0.product = p0;
						o0.setLoaded (true);
						o0.autoRebound ();
						ORF o1 = new ORF (orf1.reco, orf1.startOriginal, orf1.stopOriginal, orf1.frame);
						o1.product = p1;
						o1.setLoaded (true);
						o1.autoRebound ();
						float nacc = o0.getAccuracy () * o1.getAccuracy ();
						if (nacc > acc && o0.length + o1.length <= lgRef) {
							if (selec0 == null) {
								selec0 = o0;
								selec1 = o1;
							}
							else {
								int frO = (o0.isFragment () ? 1 : 0) + (o1.isFragment () ? 1 : 0);
								int frS = (o0.isFragment () ? 1 : 0) + (o1.isFragment () ? 1 : 0);
								int gnO = (o0.isGene () ? 1 : 0) + (o1.isGene () ? 1 : 0);
								int gnS = (selec0.isGene () ? 1 : 0) + (selec1.isGene () ? 1 : 0);
								int unO = (o0.isUncharacterized () ? 1 : 0) + (o1.isUncharacterized () ? 1 : 0);
								int unS = (o0.isUncharacterized () ? 1 : 0) + (o1.isUncharacterized () ? 1 : 0);
								if (frS > frO) {
									selec0 = o0;
									selec1 = o1;
								}
								else if (gnS < gnO) {
									selec0 = o0;
									selec1 = o1;
								}
								else if (unS > unO) {
									selec0 = o0;
									selec1 = o1;
								}
								else if (selec0.length + selec1.length < o0.length + o1.length) {
									selec0 = o0;
									selec1 = o1;
								}
							}
						}
					}
				}
				if (selec0 != null) {
					orf0.product = selec0.product;
					orf0.autoRebound ();
					orf1.product = selec1.product;
					orf1.autoRebound ();
				}
			}
			
		} catch (IOException e) {
			e.printStackTrace ();
		}
	}

}
