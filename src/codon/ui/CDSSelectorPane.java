package codon.ui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import codon.Starter;
import codon.data.CodonConfiguration;
import codon.data.ORF;
import codon.data.ProductProtein;
import codon.recognizer.OrfRecognizer;

public class CDSSelectorPane extends JPanel {
	private static final long serialVersionUID = 1L;
	
	public final static int ITEM_HEIGH = 55;

	List <ProductProtein> products = new ArrayList <ProductProtein> ();
	ORF orfRef;
	OrfRecognizer reco;
	
	int mouseY = -1;
	
	ArrayList <CDSSelectorPaneListener> listeners = new ArrayList <CDSSelectorPaneListener> ();
	public void addCDSSelectorPaneListener (CDSSelectorPaneListener l) {
		listeners.add (l);
	}
	public void removeCDSSelectorPaneListener (CDSSelectorPaneListener l) {
		listeners.remove (l);
	}
	public void fireProductChanged () {
		for (CDSSelectorPaneListener l : listeners) l.productChanged ();
	}
	
	public CDSSelectorPane (OrfRecognizer reco, ORF orf) {
		this.orfRef = orf;
		this.reco = reco;
		try {
			File f = new File (CodonConfiguration.workspace + reco.dataDestDirS + "/orfs/" + orf.id + (orf.hasSwissUpdate () ? ".swiss" : ".unip"));
			if (f.exists ()) {
				products = ProductProtein.loadAll (f);
			}
		} catch (IOException e) {
			e.printStackTrace ();
		}
		addMouseMotionListener (new MouseAdapter () {
			public void mouseMoved (MouseEvent ev) {
				mouseY = ev.getY ();
				repaint ();
			}
		});
		addMouseListener (new MouseAdapter () {
			public void mouseClicked (MouseEvent ev) {
				if (ev.getClickCount () > 1) {
					mouseY = ev.getY ();
					for (int i = 0; i < products.size (); i ++) {
						int ddy = i * ITEM_HEIGH;
						if (mouseY > ddy && mouseY < ddy + ITEM_HEIGH) {
							orfRef.product = products.get (i);
							fireProductChanged ();
							break;
						}
					}
				}
			}
		});
		setPreferredSize (new Dimension (420, products.size () * 40));
	}
	
	public void paintComponent (Graphics gg) {
		super.paintComponent (gg);
		Graphics2D g = (Graphics2D) gg;
		int ddx = 5;
		for (int i = 0; i < products.size (); i ++) {
			ProductProtein p = products.get (i);
			
			int ddy = i * ITEM_HEIGH;
			
			int boxW = 410;
			if (mouseY > ddy && mouseY < ddy + ITEM_HEIGH) {
				g.setColor (Color.LIGHT_GRAY);
				g.fillRect (0, ddy, boxW + 10, ITEM_HEIGH);
			}
			
			float [] boxes = orfRef.getMatchBox (p);
			g.setColor (Color.BLACK);
			g.drawString (p.organism, 5, ddy + 35);
			g.drawString (p.identity + "", 5, ddy + 47);
			g.drawString (p.productId + "   -  " + p.productName, 45, ddy + 47);
			
			g.fillRect ((int) (ddx + boxW * boxes [0]), ddy + 6, (int) (boxW * (boxes [1] - boxes [0])), 5);
			g.drawRect ((int) (ddx + boxW * boxes [4]), ddy + 3, (int) (boxW * (boxes [5] - boxes [4])), 17);
			
			float id = p.identity;
			if (id >= 50) g.setColor (new Color ((int) (255 * id / 100), (int) (255 * (100 - id) / 50), 0));
			else g.setColor (new Color (0, (int) (255 * id / 50), (int) (255 * (100 - id) / 100)));
			g.fillRect ((int) (ddx + boxW * boxes [2]), ddy + 13, (int) (boxW * (boxes [3] - boxes [2])), 5);
		}
	}
	
	public static void main (String [] args) {
		JFrame frame = new JFrame ();
		OrfRecognizer reco = new OrfRecognizer ();
		reco.loadFromCodonProject (new File ("workspace/ERR2348864/auto.codon"));
		ArrayList <ORF> f = reco.orfsPorFrame [0];
		ORF orfS = null;
		for (ORF orf: f) {
			if (orf.id.equals ("F0-000132-10517-10649")) {
				orfS = orf;
				orfS.start -= 18;
				orf.length += 18;
				break;
			}
		}
		JScrollPane scrollP = new JScrollPane (new CDSSelectorPane (reco, orfS));
		frame.getContentPane ().add (scrollP);
		frame.pack ();
		frame.setDefaultCloseOperation (JFrame.EXIT_ON_CLOSE);
		frame.setVisible (true);
	}
}
