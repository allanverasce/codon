package codon.ui;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.util.ArrayList;

import javax.swing.JPanel;

import codon.data.OrganismFamilly;

public class DiagramOrganisms extends JPanel {
	private static final long serialVersionUID = 1L;
	
	int max = 0;
	
	class BarOrganisms {
		private static final int BAR_WIDTH = 100;
		private static final int MARGIN_TOP = 15;
		int genes = 0;
		int proteins = 0;
		int uncharacterized = 0;
		int sum = 0;
		boolean familly = false;
		String name = "";
		public BarOrganisms (String name, boolean familly, int genes, int proteins, int uncharacterized) {
			this.genes = genes;
			this.proteins = proteins;
			this.uncharacterized = uncharacterized;
			this.familly = familly;
			this.name = name;
			sum = genes + proteins + uncharacterized;
		}
		public void paint (Graphics2D g, int pos) {
			int hRef = getHeight () - 20 - MARGIN_TOP;
//			int h = (int) (hRef * sum / (float) max);
			int hG = (int) (hRef * genes / (float) max);
			int hP = (int) (hRef * proteins / (float) max);
			int hU = (int) (hRef * uncharacterized / (float) max);
			g.setColor (Color.BLACK);
			g.drawString (name, pos * BAR_WIDTH, MARGIN_TOP + hRef + 15);
			g.fillRect (pos * BAR_WIDTH, MARGIN_TOP + hRef - hG, BAR_WIDTH, hG);
			g.setColor (Color.DARK_GRAY);
			g.fillRect (pos * BAR_WIDTH, MARGIN_TOP + hRef - hG - hP, BAR_WIDTH, hP);
			g.fillRect (pos * BAR_WIDTH, MARGIN_TOP + hRef - hG - hP - hU, BAR_WIDTH, hU);
		}
	}
	
	ArrayList <OrganismFamilly> famillies;
	ArrayList <BarOrganisms> bars = new ArrayList <DiagramOrganisms.BarOrganisms> ();
	
	public DiagramOrganisms (ArrayList <OrganismFamilly> famillies) {
		setLayout (null);
		this.famillies = famillies;
		for (int i = 0; i < famillies.size (); i++) {
			OrganismFamilly fam = famillies.get (i);
			bars.add (new BarOrganisms (fam.famillyName, true, fam.genes, fam.proteins, fam.uncharacterized));
		}
	}
	
	protected void paintComponent (Graphics arg0) {
		Graphics2D g = (Graphics2D) arg0;
		int pos = 0;
		for (BarOrganisms bar: bars) {
			bar.paint (g, pos);
			pos ++;
		}
	}
}
