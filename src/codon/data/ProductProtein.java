package codon.data;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class ProductProtein {
	public String productName = "";
	public String productId = "";
	public String gene = "";
	public String entryId = "no entry";
	public float identity = 0;
	public float score = 0;
	public boolean blasted = false;
	public boolean hasEntries = false;
	public boolean fragment = false;
	public boolean swiss = false;
	public boolean uncharacterized = false; 
	
	private int startMatchSeq = 0;
	private int endMatchSeq = 0;
	private int startQuerySeq = 0;
	private int endQuerySeq = 0;
	public String matchingSequence = "";
	
	public String goTerm = "";
	public String kegg = "";
	
	public String comment = "";
	public String organism = "";
	
	public boolean hasMatchingSequence () {
		return matchingSequence.length () > 0;
	}
	
	public ProductProtein () {}
	
	public ProductProtein (String dataFileS, String id, boolean fromSwiss) {
		File f = new File (CodonConfiguration.workspace + dataFileS + "/orfs/" + id + (fromSwiss ? ".swiss" : ".unip"));
		swiss = fromSwiss;
		if (!f.exists ()) return;
		blasted = true;
		try {
			BufferedReader reader = new BufferedReader (new FileReader (f));
			loadProduct (reader);
			reader.close ();
		} catch (FileNotFoundException e) {
			e.printStackTrace ();
		} catch (IOException e) {
			System.err.println ("Error while reading " + f.getAbsolutePath ());
			e.printStackTrace ();
		}
	}
	
	public static List <ProductProtein> loadAll (String f) throws IOException {
		return loadAll (new File (f));
	}
	
	public static List <ProductProtein> loadAll (File f) throws IOException {
		ArrayList<ProductProtein> products = new ArrayList <ProductProtein> ();
		BufferedReader reader = new BufferedReader (new FileReader (f));
		ProductProtein product = new ProductProtein ();
		product.loadProduct (reader);
		while (product.hasEntries) {
			products.add (product);
			product = new ProductProtein ();
			product.loadProduct (reader);
		};
		reader.close ();
		return products;
	}
	
	public static void save (ProductProtein product, String f) throws IOException {
		ArrayList <ProductProtein> pds = new ArrayList <ProductProtein> ();
		pds.add (product);
		saveAll (pds, f);
	}
	
	public static void saveAll (List <ProductProtein> products, String f) throws IOException {
		saveAll (products, new File (f));
	}
	
	public static void saveAll (List <ProductProtein> products, File f) throws IOException {
		File dir = f.getParentFile ();
		if (!dir.exists ()) dir.mkdirs ();
		BufferedWriter writer = new BufferedWriter (new FileWriter (f));
		for (ProductProtein p: products) {
			p.save (writer);
		}		
		writer.close ();
	}
	
	public ProductProtein copy () {
		ProductProtein copy = new ProductProtein ();
		copy.productName = productName;
		copy.productId = productId;
		copy.gene = gene;
		copy.entryId = entryId;
		copy.identity = identity;
		copy.score = score;
		copy.blasted = blasted;
		copy.hasEntries = hasEntries;
		copy.fragment = fragment;
		copy.swiss = swiss;
		copy.uncharacterized = uncharacterized; 
		copy.startMatchSeq = startMatchSeq;
		copy.endMatchSeq = endMatchSeq;
		copy.startQuerySeq = startQuerySeq;
		copy.endQuerySeq = endQuerySeq;
		copy.matchingSequence = matchingSequence;
		copy.goTerm = goTerm;
		copy.kegg = kegg;
		copy.comment = comment;
		copy.organism = organism;
		return copy;
	}
	
	public void loadProduct (BufferedReader reader) throws IOException {
		blasted = true;
		entryId = reader.readLine ();
		if (entryId != null) {
			hasEntries = true;
		}
		else {
			return;
		}
		productName = reader.readLine ();
		if (productName == null) throw new IOException ();
		fragment = productName.contains ("Fragment");
		productId = reader.readLine ();
		gene = reader.readLine ();
		uncharacterized = productName.contains ("Uncharacterized") && !isGene ();
		identity = Float.parseFloat (reader.readLine ());
		score = Float.parseFloat (reader.readLine ());
		setStartMatchSeq (Integer.parseInt (reader.readLine ()));
		setEndMatchSeq (Integer.parseInt (reader.readLine ()));
		setStartQuerySeq (Integer.parseInt (reader.readLine ()));
		setEndQuerySeq (Integer.parseInt (reader.readLine ()));
		matchingSequence = reader.readLine ();
		goTerm = reader.readLine ();
		kegg = reader.readLine ();
		organism = reader.readLine ();
		organism = organism == null ? "" : organism;
	}
	
	public void save (String dataFileS, String id, boolean fromSwiss) {
		try {
			File f = new File (CodonConfiguration.workspace + dataFileS + "/orfs/" + id + (fromSwiss ? ".swiss" : ".unip"));
			BufferedWriter writer = new BufferedWriter (new FileWriter (f));
			writer.write (entryId); writer.newLine ();
			writer.write (productName); writer.newLine ();
			writer.write (productId); writer.newLine ();
			writer.write (gene); writer.newLine ();
			writer.write (identity + ""); writer.newLine ();
			writer.write (score + ""); writer.newLine ();
			writer.write ((getStartMatchSeq ()) + ""); writer.newLine ();
			writer.write (getEndMatchSeq () + ""); writer.newLine ();
			writer.write ((getStartQuerySeq ()) + ""); writer.newLine ();
			writer.write ((getEndQuerySeq ()) + ""); writer.newLine ();
			writer.write (matchingSequence); writer.newLine ();
			writer.write (goTerm); writer.newLine ();
			writer.write (kegg); writer.newLine ();
			writer.write (organism); writer.newLine ();
			writer.close ();
		} catch (IOException e) {
			e.printStackTrace ();
		}
	}
	
	public void save  (BufferedWriter writer) throws IOException {
		writer.write (entryId); writer.newLine ();
		writer.write (productName); writer.newLine ();
		writer.write (productId); writer.newLine ();
		writer.write (gene); writer.newLine ();
		writer.write (identity + ""); writer.newLine ();
		writer.write (score + ""); writer.newLine ();
		writer.write ((getStartMatchSeq ()) + ""); writer.newLine ();
		writer.write (getEndMatchSeq () + ""); writer.newLine ();
		writer.write ((getStartQuerySeq ()) + ""); writer.newLine ();
		writer.write ((getEndQuerySeq ()) + ""); writer.newLine ();
		writer.write (matchingSequence); writer.newLine ();
		writer.write (goTerm); writer.newLine ();
		writer.write (kegg); writer.newLine ();
		writer.write (organism); writer.newLine ();
	}
	
	public boolean isGene () {
		return gene != null && gene.length () > 0;
	}
	
	public float getSpecificity () {
		return identity * (getEndMatchSeq () - getStartMatchSeq ()) / matchingSequence.length ();
	}
	
	public boolean isUncharacterized () {
		return uncharacterized;
	}
	
	public String toString () {
		return 	productName + "�" + productId + "�" + gene + "�" + entryId + "�" + identity + "�" + score + "�" + blasted + "�" + hasEntries + "�" 
				+ fragment + "�" + swiss + "�" + uncharacterized + "�" + getStartMatchSeq () + "�" + getEndMatchSeq () + "�" + getStartQuerySeq () + "�" + getEndQuerySeq () + "�" 
				+ matchingSequence + "�" + goTerm + "�" + kegg + "�" + organism + "�" + comment;
	}

	public static ProductProtein createFromStrings (String [] args) {
		ProductProtein product = new ProductProtein ();
		product.productName = args [8];
		product.productId = args [9];
		product.gene = args [10];
		product.entryId = args [11];
		product.identity = Float.parseFloat (args [12]);
		product.score = Float.parseFloat (args [13]);
		product.blasted = Boolean.parseBoolean (args [14]);
		product.hasEntries = Boolean.parseBoolean (args [15]);
		product.fragment = Boolean.parseBoolean (args [16]);
		product.swiss = Boolean.parseBoolean (args [17]);
		product.uncharacterized = Boolean.parseBoolean (args [18]);
		product.setStartMatchSeq(Integer.parseInt (args [19]));
		product.setEndMatchSeq(Integer.parseInt (args [20]));
		product.setStartQuerySeq(Integer.parseInt (args [21]));
		product.setEndQuerySeq(Integer.parseInt (args [22]));
		if (args.length > 23) product.matchingSequence = args [23];
		if (args.length > 24) product.goTerm = args [24];
		if (args.length > 25) product.kegg = args [25];
		if (args.length > 26) product.organism  = args [26];
		if (args.length > 27) product.comment = args [27];
		return product;
	}

	public int getStartQuerySeq () {
		return startQuerySeq;
	}

	public void setStartQuerySeq (int startQuerySeq) {
		this.startQuerySeq = startQuerySeq;
	}

	public int getStartMatchSeq () {
		return startMatchSeq;
	}

	public void setStartMatchSeq (int startMatchSeq) {
		this.startMatchSeq = startMatchSeq;
	}

	public int getEndQuerySeq () {
		return endQuerySeq;
	}

	public void setEndQuerySeq (int endQuerySeq) {
		this.endQuerySeq = endQuerySeq;
	}

	public int getEndMatchSeq () {
		return endMatchSeq;
	}

	public void setEndMatchSeq (int endMatchSeq) {
		this.endMatchSeq = endMatchSeq;
	}
	
}
