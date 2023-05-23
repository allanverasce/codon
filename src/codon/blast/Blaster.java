package codon.blast;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;

import codon.data.CodonConfiguration;
import codon.data.ORF;
import codon.data.ProductProtein;
import codon.recognizer.OrfRecognizer;
import uk.ac.ebi.kraken.interfaces.uniprot.DatabaseCrossReference;
import uk.ac.ebi.kraken.interfaces.uniprot.Gene;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;
import uk.ac.ebi.kraken.interfaces.uniprot.dbx.go.Go;
import uk.ac.ebi.uniprot.dataservice.client.Client;
import uk.ac.ebi.uniprot.dataservice.client.ServiceFactory;
import uk.ac.ebi.uniprot.dataservice.client.alignment.blast.Alignment;
import uk.ac.ebi.uniprot.dataservice.client.alignment.blast.BlastInput;
import uk.ac.ebi.uniprot.dataservice.client.alignment.blast.BlastResult;
import uk.ac.ebi.uniprot.dataservice.client.alignment.blast.UniProtBlastService;
import uk.ac.ebi.uniprot.dataservice.client.alignment.blast.UniProtHit;
import uk.ac.ebi.uniprot.dataservice.client.alignment.blast.input.DatabaseOption;
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtService;

public class Blaster {
	
	public boolean receiving = false;
	public boolean stopped = false;
	
	public static int ongoingBlast = 0;
	public static synchronized void enteringBlast () {
		ongoingBlast ++;
	}
	public static synchronized void leavingBlast () {
		ongoingBlast --;
	}
	
	public static int ongoingEntry = 0;
	public static synchronized void enteringEntry () {
		ongoingEntry ++;
	}
	public static synchronized void leavingEntry () {
		ongoingEntry --;
	}
		
	public static LinkedList <Boolean> hitsSwiss = new LinkedList <> ();
	public static LinkedList <Long> hitsTimeSwiss = new LinkedList <> ();
	public static LinkedList <Long> hitsTimeUniprot = new LinkedList <> ();
	public static int averageWindowSize = 60;
	public static boolean acceleratedWithSwissProt = true;
	
	public static long averageSwissHitTime () {
		if (hitsTimeSwiss.size () == 0) return 0;
		long averageTimeSwiss = 0;
		for (Long tm: hitsTimeSwiss) averageTimeSwiss += tm;
		averageTimeSwiss /= hitsTimeSwiss.size ();
		return averageTimeSwiss;
	}
	
	public static long averageUniprotHitTime () {
		long averageTimeUniprot = 0;
		if (hitsTimeUniprot.size () == 0) return 0;
		for (Long tm: hitsTimeUniprot) averageTimeUniprot += tm;
		averageTimeUniprot /= hitsTimeUniprot.size ();
		return averageTimeUniprot;
	}
	
	public static int successedHitCount () {
		int successedHitCount = 0;
		for (boolean hit: hitsSwiss) {
			if (hit) successedHitCount ++;
		}
		return successedHitCount;
	}
	
	public float blast (ORF orf, String baseDir) {
		if (baseDir != null) {
			if (new File (CodonConfiguration.workspace + baseDir + "/orfs/" + orf.id + ".unip").exists ()) return 100;
			if (new File (CodonConfiguration.workspace + baseDir + "/orfs/" + orf.id + ".swiss").exists ()) return 100;
			File dir = new File (CodonConfiguration.workspace + baseDir + "/orfs/");
			if (!dir.exists ()) dir.mkdirs ();
		}
		this.baseDir = baseDir;
//		return blast (orf, DatabaseOption.SWISSPROT);
		if (acceleratedWithSwissProt) {
			int successedHitsCount = successedHitCount ();
			long averageTimeSwiss = averageSwissHitTime ();
			long averageTimeUniprot = averageUniprotHitTime ();
			
			long tm1 = successedHitsCount * averageTimeSwiss;
			long tm2 = (averageWindowSize - successedHitsCount) * (averageTimeSwiss + averageTimeUniprot);
			long tm3 = averageWindowSize * averageTimeUniprot;
					
			if (hitsTimeUniprot.size () < averageWindowSize / 2 || hitsTimeSwiss.size () < averageWindowSize || hitsTimeSwiss.size () == averageWindowSize && tm1 + tm2 <= tm3) {
				long tm = System.currentTimeMillis ();
				swissOngoing = true;
				blast (orf, DatabaseOption.SWISSPROT);
				swissOngoing = false;
				if (maxAcc != -1) {
					if (hitsSwiss.size () == averageWindowSize) {
						hitsTimeSwiss.remove (0);
						hitsSwiss.remove (0);
					}
					hitsSwiss.add (maxAcc >= 90);
					hitsTimeSwiss.add (System.currentTimeMillis () - tm);
				}
				if (maxAcc > 90) {
					swiss = true;
				}
			}
			else {
				acceleratedWithSwissProt = false;
			}
		}
		if (maxAcc < 90) {
			long tm = System.currentTimeMillis ();
			blast (orf, DatabaseOption.UNIPROTKB);
			if (maxAcc != -1) {
				if (hitsTimeUniprot.size () == averageWindowSize) hitsTimeUniprot.remove (0);
				hitsTimeUniprot.add (System.currentTimeMillis () - tm);
			}
		}
		if (maxAcc != -1) save ();
		return maxAcc;
	}

	ArrayList <ProductProtein> products = new ArrayList <> ();
	
	String baseDir = null;
	String id = null;
	float maxAcc = -1;
	boolean swiss = false;
	
	boolean swissOngoing = false;
	
	public float blast (ORF orf, DatabaseOption database) {
		String sequence = orf.getOriginalAmidoSequence ();
		sequence = sequence.replaceAll ("\\+", "");
		sequence = sequence.replaceAll ("\\*", "");
		sequence = sequence.replaceAll ("\\#", "");
//		if (sequence.length () * 3 < ORF.MIN_LENGHT) {
//			System.err.println (orf.id + ":" + sequence);
//			return true;
//		}
		receiving = false;
		this.id = orf.id;
		
		String querySequence = orf.id + "\n" + sequence;
			
		ServiceFactory serviceFactoryInstance = Client.getServiceFactoryInstance ();
	    UniProtBlastService uniProtBlastService = serviceFactoryInstance.getUniProtBlastService ();
	    enteringBlast ();
	    uniProtBlastService.start ();
	    BlastInput input = new BlastInput.Builder (database, querySequence).build ();
	    CompletableFuture <BlastResult <UniProtHit>> resultFuture = uniProtBlastService.runBlast (input);

	    try {
	    	products.addAll (retrieveResults (resultFuture));
	    	if (stopped) {
	    		uniProtBlastService.stop ();
	    		return -1;
	    	}
//            if (products.size () > 0) {
	    	float max = 0;
			for (ProductProtein p: products) {
				ORF orf1 = new ORF (orf.reco, orf.startOriginal, orf.stopOriginal, orf.frame);
				orf1.product = p;
				orf1.setLoaded (true);
				orf1.autoRebound ();
				float acc = orf1.getAccuracy ();
				if (max < acc) max = acc;
			}
			
			maxAcc = max;
        	try {
				Collections.sort (products, new Comparator <ProductProtein> () {
					public int compare (ProductProtein p1, ProductProtein p2) {
						ORF orf1 = new ORF (orf.reco, orf.startOriginal, orf.stopOriginal, orf.frame);
						ORF orf2 = new ORF (orf.reco, orf.startOriginal, orf.stopOriginal, orf.frame);
						orf1.product = p1;
						orf1.setLoaded (true);
						orf1.autoRebound ();
						orf2.product = p2;
						orf2.setLoaded (true);
						orf2.autoRebound ();
						
						float acc1 = orf1.getAccuracy ();
						float acc2 = orf2.getAccuracy ();
						
//							if (acc1 < 20 && acc2 < 20) return acc1 > acc2 ? -1 : (acc1 == acc2 ? 0 : 1);
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
        	}
        	catch (Exception e) {
        		System.err.println ("Error while sorting product for " + orf.id);
			}
        } catch (ExecutionException e) {
        	e.printStackTrace ();
        	maxAcc = -1;
        } catch (InterruptedException e) {
        	e.printStackTrace ();
        	maxAcc = -1;
		} 
	    uniProtBlastService.stop ();
	    leavingBlast ();
	    return maxAcc;
	}
	
	public void save () {
		if (baseDir != null) {
			try {
				ProductProtein.saveAll (products, CodonConfiguration.workspace + baseDir + "/orfs/" + id + (swiss ? ".swiss" : ".unip"));
			} catch (IOException e) {
				e.printStackTrace ();
			}
		}
	}
	
	private UniProtEntry retrieveEntry (String entryId) {
        ServiceFactory serviceFactoryInstance = Client.getServiceFactoryInstance ();
        UniProtService uniProtService = serviceFactoryInstance.getUniProtQueryService ();
        UniProtEntry entry = null;
        try {
            uniProtService.start ();
            enteringEntry ();
            entry = uniProtService.getEntry (entryId);
        } catch (Exception e) {
        	e.printStackTrace ();
        } 
        leavingEntry ();
        return entry;
    }
	
	private ArrayList <ProductProtein> retrieveResults (CompletableFuture <BlastResult <UniProtHit>> resultFuture) throws InterruptedException, ExecutionException {
		BlastResult <UniProtHit> blastResult = resultFuture.get ();
        ArrayList <ProductProtein> products = new ArrayList <> (); 
        for (UniProtHit hit : blastResult.hits ()) {
        	if (stopped) return null;
        	ProductProtein product = new ProductProtein ();
        	receiving = true;
        	if (hit == null) continue;
        	
        	UniProtEntry entry = hit.getEntry ();
        	
        	if (hit.getSummary () != null && entry == null) {
        		String entryId = hit.getSummary ().getEntryId ();
        		if (entryId.contains ("_")) {
        			entryId = entryId.substring (0, entryId.indexOf ('_'));
        			entry = retrieveEntry (entryId);
        		}
        	}
			if (entry == null || hit.getEntry ().getPrimaryUniProtAccession () == null || hit.getSummary () == null) continue;
			
			String desc = hit.getSummary ().getDescription ();
			
			product.entryId = entry.getPrimaryUniProtAccession ().getValue ();
			product.productName = desc.substring (0, desc.indexOf (" OS="));
			product.uncharacterized = product.productName.contains ("Uncharacterized");
			product.fragment = product.productName.contains ("Fragment");
			
			String productId = "";
			if (desc.contains ("GN=")) {
				if (desc.contains (" PE=")) productId = desc.substring (desc.indexOf ("GN=") + 3, desc.indexOf (" PE="));
				else productId = desc.substring (desc.indexOf ("GN=") + 3);
			}
			else productId = "";
			product.productId = productId;
			
			List <Gene> genes = entry.getGenes ();
			String gene = "";
        	if (genes.size () > 0) gene = genes.get (0).getGeneName ().getValue ();
        	product.gene = gene;
			
			List <Alignment> aligns =  hit.getSummary ().getAlignments ();
			
			for (Alignment al : aligns) {
				product.identity = al.getIdentity ();
				product.score = al.getScore ();
				product.setStartMatchSeq (al.getStartMatchSeq () - 1);
				product.setEndMatchSeq (al.getEndMatchSeq ());
				product.setStartQuerySeq (al.getStartQuerySeq () - 2);
				product.setEndQuerySeq (al.getEndQuerySeq () - 1);
				break;
			}
			
			String matchingSequence = entry.getSequence ().getValue ();
			
			product.matchingSequence = matchingSequence;
			List <Go> gos = entry.getGoTerms ();
            for (Go go: gos) {
            	product.goTerm += go.getGoId () + ":" + go.getDescription () + ": " + go.getGoTerm () + "|";
            }
        	Collection <DatabaseCrossReference> dts = entry.getDatabaseCrossReferences ();
            for (DatabaseCrossReference dt: dts) {
            	if (dt.getDatabase ().getName ().contains ("KEGG")) {
            		product.kegg += dt.getPrimaryId ().getValue () + "|";
             	}
            }
            product.organism = entry.getOrganism ().toString ();
            products.add (product);
        }
        return products;
	}
	
	public static void reSort (OrfRecognizer reco, ORF orf) throws IOException {
		orf.loadProduct ();
		if (orf.isBlasted ()) {
			List <ProductProtein> products = ProductProtein.loadAll (CodonConfiguration.workspace + reco.dataDestDirS + "/orfs/" + orf.id + ".unip");
			float max = 0;
			for (ProductProtein p: products) {
				ORF orf1 = new ORF (orf.reco, orf.startOriginal, orf.stopOriginal, orf.frame);
				orf1.product = p;
				orf1.setLoaded (true);
				orf1.autoRebound ();
				float acc = orf1.getAccuracy ();
				if (max < acc) max = acc;
			}
			final float maxAcc = max;
			Collections.sort (products, new Comparator <ProductProtein> () {
				public int compare (ProductProtein p1, ProductProtein p2) {
					ORF orf1 = new ORF (orf.reco, orf.startOriginal, orf.stopOriginal, orf.frame);
					ORF orf2 = new ORF (orf.reco, orf.startOriginal, orf.stopOriginal, orf.frame);
					orf1.product = p1;
					orf1.setLoaded (true);
					orf1.autoRebound ();
					orf2.product = p2;
					orf2.setLoaded (true);
					orf2.autoRebound ();
					
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
					if (orf1.length < orf2.length) return -1;
					
					if (orf1.length > orf2.length) return 1;
					if (acc1 == acc2) return 0;
					if (acc1 > acc2) return -1;
					if (acc1 < acc2) return 1;
					return 0;
				}
			});
			ProductProtein.saveAll (products, CodonConfiguration.workspace + reco.dataDestDirS + "/orfs/" + orf.id + ".unip");
		}
		if (orf.hasSwissUpdate ()) {
			List <ProductProtein> products = ProductProtein.loadAll (CodonConfiguration.workspace + reco.dataDestDirS + "/orfs/" + orf.id + ".swiss");
			float max = 0;
			for (ProductProtein p: products) {
				ORF orf1 = new ORF (orf.reco, orf.startOriginal, orf.stopOriginal, orf.frame);
				orf1.product = p;
				orf1.setLoaded (true);
				orf1.autoRebound ();
				float acc = orf1.getAccuracy ();
				if (max < acc) max = acc;
			}
			final float maxAcc = max;
			Collections.sort (products, new Comparator <ProductProtein> () {
				public int compare (ProductProtein p1, ProductProtein p2) {
					ORF orf1 = new ORF (orf.reco, orf.startOriginal, orf.stopOriginal, orf.frame);
					ORF orf2 = new ORF (orf.reco, orf.startOriginal, orf.stopOriginal, orf.frame);
					orf1.product = p1;
					orf1.setLoaded (true);
					orf1.autoRebound ();
					orf2.product = p2;
					orf2.setLoaded (true);
					orf2.autoRebound ();
					
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
					if (orf1.length < orf2.length) return -1;
					
					if (orf1.length > orf2.length) return 1;
					if (acc1 == acc2) return 0;
					if (acc1 > acc2) return -1;
					if (acc1 < acc2) return 1;
					return 0;
				}
			});
			ProductProtein.saveAll (products, CodonConfiguration.workspace + reco.dataDestDirS + "/orfs/" + orf.id + ".swiss");
		}
	}
	
	public static void reSort (OrfRecognizer reco) throws IOException {
//		int k = 0;
		
		File fls [] = new File (CodonConfiguration.workspace + reco.dataDestDirS + "/orfs/").listFiles ();
		for (int i = 0; i < fls.length; i++) {
			if (fls [i].getName ().contains (".unip")) {
				String n = fls [i].getName ().substring (1);
				int f = Integer.parseInt (n.substring (0, 1));
				n = n.substring (n.indexOf ("-", 2) + 1);
				int start = Integer.parseInt (n.substring (0, n.indexOf ("-")));
				n = n.substring (n.indexOf ("-") + 1);
				int stop = Integer.parseInt (n.substring (0, n.indexOf (".")));
				ORF orf = new ORF (reco, start, stop, f);
				if (orf.id.equals ("F0-001542-482-2024")) {
//					System.err.println (++k + " " + orf.id );
					try {
						reSort (reco, orf);
					}
					catch (Exception e) {
						e.printStackTrace ();
						System.out.println ("Error in " + orf.id);
					}
				}
			}
		}
	}
	
//	public static void main (String [] args) {
//		OrfRecognizer reco = new OrfRecognizer ();
//		reco.load ("workspace/PALESB58.fasta");
//		reco.dataDestFileS = "PALESB58";
//		try {
//			reSort (reco);
//		} catch (IOException e) {
//			e.printStackTrace ();
//		}
//	}
}
