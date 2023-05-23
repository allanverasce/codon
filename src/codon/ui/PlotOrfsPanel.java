package codon.ui;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.Stroke;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.font.TextAttribute;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Vector;

import javax.swing.JPanel;
import javax.swing.SwingUtilities;

import codon.data.Codon;
import codon.data.CodonConfiguration;
import codon.data.GeneFeature;
import codon.data.ORF;
import codon.recognizer.OrfRecognizer;

import java.awt.Font;

public class PlotOrfsPanel extends JPanel {
	private static final long serialVersionUID = 1L;

//	public static final Color COR_NORMAL_WAY = Color.BLUE;
	public static final Color COR_NORMAL_WAY = new Color (64, 64, 64);
	public static final Color COR_NORMAL_WAY_REMOVED = new Color (64, 64, 64, 20);
	public static final Color COR_REVERSE_WAY_REMOVED = new Color (64, 64, 64, 20);
//	public static final Color COR_REVERSE_WAY = Color.RED;
	public static final Color COR_REVERSE_WAY = new Color (64, 64, 64);
	
	public static final Color FORWARD_BACKGROUND_COLOR = new Color (205, 210, 205);
	public static final Color REVERSE_BACKGROUND_COLOR = new Color (210, 205, 205);
	
	public static final Color FRAME_SEPARATION_COLOR = new Color (235, 235, 235);
	public static final Color LINE_SEPARATION_COLOR = new Color (215, 215, 215);
	public static final Color LINE_SEPARATION_TOP_BORDER_COLOR = Color.WHITE;
	public static final Color LINE_SEPARATION_BOTTOM_BORDER_COLOR = Color.DARK_GRAY;
	public static final Color LINE_SEPARATION_TEXT_COLOR = Color.BLACK;
	
	public static final Color INTERGENE_AREA_COLOR = new Color (245, 245, 245);
	public static final Color OVERLAP_AREA_COLOR = new Color (170, 170, 170);
	public static final Color RRNA_AREA_COLOR = new Color (0, 0, 128);
	public static final Color TRNA_AREA_COLOR = new Color (50, 50, 128);
	
	
	public static final Color COR_STOP_CODON = new Color (0, 0, 180);
	public static final Color COR_START_CODON = new Color (0, 180, 0);
    public static final Color COR_VALINA = Color.YELLOW;
	public static final Color COR_SELECT = new Color (255, 255, 255, 230);
	public static final Color COR_BCK_SELECT = new Color (230, 230, 230, 255);
	
	public static final Stroke SELECT_STROKE = new BasicStroke (7);
	public static final Stroke NORMAL_STROKE = new BasicStroke (1);
	final static BasicStroke DASHED_STROKE = new BasicStroke (1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, new float [] {10.0f}, 0.0f);
	
	OrfRecognizer reco;
	
	int scaleY = 14;
	int offsetY = 0;
	
	SequenceEditor sequenceEditor = null;
	
	class ORFBox extends Rectangle {
		private static final long serialVersionUID = 1L;
		ORF orf;
		public ORFBox (ORF par, int x, int y, int w, int h) {
			super (x, y, w, h);
			this.orf = par;
		}
	}
	ArrayList <ORFBox> orfBoxes = new ArrayList <> ();
	private ORFBox selected = null;
	boolean computeBox = true;
	PlotOrfsFrame plotPairsJFrame;
	
	public PlotOrfsPanel (OrfRecognizer reco, PlotOrfsFrame plotPairsJFrame) {
		this.reco = reco;
		this.plotPairsJFrame = plotPairsJFrame;
		setFocusable (true);
		grabFocus ();
		addMouseWheelListener (new MouseWheelListener () {
			public void mouseWheelMoved (MouseWheelEvent e) {
				setOffsetY (offsetY + e.getWheelRotation () * getHeight () / ((int) (scaleY * 7.5 * 2)), true);
			}
		});
		addMouseListener (new MouseAdapter () {
			public void mousePressed (MouseEvent ev) {
				grabFocus ();
//				((JComponent) plotPairsJFrame.getContentPane ()).grabFocus ();
				if (ev.getClickCount () == 2) {
					if (getSelected () == null) return;
					getSelected ().orf.loadProduct ();
					new CDSViewDialog (PlotOrfsPanel.this, getSelected ().orf, plotPairsJFrame.menuViewAcidoSequence.isSelected ());
				}
				else {
					for (int i = 0; i < orfBoxes.size (); i++) {
						if (orfBoxes.get (i).contains (ev.getX (), ev.getY ())) {
							if (orfBoxes.get (i) != getSelected ()) {
								setSelected (orfBoxes.get (i));
								repaint ();
							}
							return;
						}
					}
					
					int space = (int) (7.5 * scaleY);
					int line = ev.getY () / space;
					int frame = (int) ((ev.getY () - space * line) / scaleY);
					int index = offsetY  * getWidth () + line * getWidth () + ev.getX ();
					
					if (frame == 0) {
						sequenceEditor.setFirstIndex (index);
						repaint ();
						sequenceEditor.repaint ();
					}
					if (ev.isControlDown ()) {
						if (frame > 0 && frame < 7) {
							reco.createOrf (frame - 1, index);
							computeBox = true;
							repaint ();
						}
					}
					else {
						if (getSelected () != null) {
							setSelected (null);
							repaint ();
						}
					}
				}
			}
		});
		addKeyListener (new KeyAdapter () {
			public void keyPressed (KeyEvent ev) {
				if (ev.getKeyCode () == KeyEvent.VK_PAGE_UP) {
					setOffsetY (offsetY - getHeight () / (scaleY * 7) - 1, true);
				}
				if (ev.getKeyCode () == KeyEvent.VK_PAGE_DOWN) {
					setOffsetY (offsetY + getHeight () / (scaleY * 7) - 1, true);
				}
				if (ev.getKeyCode () == KeyEvent.VK_UP) {
					setOffsetY (offsetY - 1, true);
				}
				if (ev.getKeyCode () == KeyEvent.VK_DOWN) {
					setOffsetY (offsetY + 1, true);
				}
//				if (ev.getKeyCode () == KeyEvent.VK_ADD) {
//					scaleY += 2;
//					repaint ();
//				}
//				if (ev.getKeyCode () == KeyEvent.VK_SUBTRACT) {
//					if (scaleY > 14) scaleY -= 2;
//					repaint ();
//				}
			}
		});
		addComponentListener (new ComponentAdapter () {
			public void componentResized (ComponentEvent e) {
				recompute ();
			}
		});
	}

	private void setOffsetY (int offset, boolean updateScrollBar) {
		this.offsetY = offset;
		if (offsetY < 0) offsetY = 0;
		if (updateScrollBar) {
			plotPairsJFrame.scroll.setValue (offsetY * getWidth ());
		}
		recompute ();
	}
	
	public int getVisibleIndex () {
		return offsetY * 7 * scaleY;
	}
	
	public void setSelectedPosition (int pos, boolean updateScrollBar) {
		setOffsetY (pos / getWidth (), updateScrollBar);
		if (sequenceEditor != null) {
			sequenceEditor.setFirstIndex (pos);
			sequenceEditor.repaint ();
		}
	}
	
	public void recompute () {
		computeBox = true;
		mediaMovelGC = null;
		setSelected (null);
		repaint ();
	}

	Vector <Float> mediaMovelGC = null;
	
	protected void paintComponent (Graphics gg) {
		super.paintComponent (gg);
		Graphics2D g = (Graphics2D) gg;
		g.setRenderingHint (RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_LCD_HRGB);
		int w = getWidth ();
		int h = getHeight ();

		int space = (int) (7.5 * scaleY);
		int line_separation_height = scaleY;
		int line_separation_offset = 0;
		
		int frames_height = 6 * scaleY;
		int frame_height = scaleY;
		int frames_offset = scaleY;
		
		g.setColor (FORWARD_BACKGROUND_COLOR);
		for (int i = 0; i < h; i += space) {
			g.fillRect (0, i + frames_offset, w, space);
		}
		
		g.setColor (REVERSE_BACKGROUND_COLOR);
		for (int i = 0; i < h; i += space) {
			g.fillRect (0, i + frames_offset + 3 * frame_height, w, frame_height * 3);
		}
		
		g.setColor (LINE_SEPARATION_COLOR);
		for (int i = 0; i < h; i += space) {
			g.fillRect (0, i + line_separation_offset, w, line_separation_height);
		}
			
		g.setColor (LINE_SEPARATION_TEXT_COLOR);
		for (int i = 0; i < h; i += space) {
			int v = (offsetY + i / space) * w;
			g.drawString (v + "", 5, line_separation_offset + frame_height + i - 2);
		}
		g.setColor (LINE_SEPARATION_TOP_BORDER_COLOR);
		for (int i = 0; i < h; i += space) {
			g.drawLine (0, i + line_separation_offset, w, i + line_separation_offset);
		}
		g.setColor (LINE_SEPARATION_BOTTOM_BORDER_COLOR);
		for (int i = 0; i < h; i += space) {
			g.drawLine (0, i + line_separation_offset + frame_height, w, i + line_separation_offset + frame_height);
		}
		
		if (plotPairsJFrame.menuViewIntergene.isSelected ()) {
			List <Codon> intergeneSpans = reco.intergenicRegions;
			List <Codon> superpositionSpans = reco.overlapSpans;
			g.setColor (INTERGENE_AREA_COLOR);
			for (int i = 0; i < intergeneSpans.size () - 1; i += 2) {
				int start = intergeneSpans.get (i).index;
				int stop = intergeneSpans.get (i + 1).index;
				int y1 = line_separation_offset + space * (start / w) - offsetY * space;
				int x1 = start % w;
				int y2 = line_separation_offset + space * (stop / w) - offsetY * space;
				if (y2 < 0) continue;
				int x2 = stop % w;
				if (y2 != y1) {
					g.fillRect (x1, y1, w - x1, line_separation_height);
					for (int j = y1 + space; j < y2; j += space) {
						g.fillRect (0, j, w, line_separation_height);
					}
					g.fillRect (0, y2, x2, line_separation_height);
				}
				else {
					g.fillRect (x1, y1, x2 - x1, line_separation_height);
				}
				if (y1 > h) break;
			}
			g.setColor (OVERLAP_AREA_COLOR);
			for (int i = 0; i < superpositionSpans.size () - 1; i += 2) {
				int start = superpositionSpans.get (i).index;
				int stop = superpositionSpans.get (i + 1).index;
				int y1 = line_separation_offset + space * (start / w) - offsetY * space;
				int x1 = start % w;
				int y2 = line_separation_offset + space * (stop / w) - offsetY * space;
				if (y2 < 0) continue;
				int x2 = stop % w;
				if (y2 != y1) {
					g.fillRect (x1, y1, w - x1, line_separation_height);
					for (int j = y1 + space; j < y2; j += space) {
						g.fillRect (0, j, w, line_separation_height);
					}
					g.fillRect (0, y2, x2, line_separation_height);
				}
				else {
					g.fillRect (x1, y1, x2 - x1, line_separation_height);
				}
				if (y1 > h) break;
			}
		}
		List <ORF> rRNA = reco.rnaOrfs;
		g.setColor (RRNA_AREA_COLOR);
		for (int i = 0; i < rRNA.size (); i ++) {
			int start = rRNA.get (i).start;
			int stop = rRNA.get (i).stop;
			int y1 = line_separation_offset + space * (start / w) - offsetY * space;
			int x1 = start % w;
			int y2 = line_separation_offset + space * (stop / w) - offsetY * space;
			if (y2 < 0) continue;
			int x2 = stop % w;
			if (y2 != y1) {
				g.fillRect (x1, y1, w - x1, line_separation_height);
				for (int j = y1 + space; j < y2; j += space) {
					g.fillRect (0, j, w, line_separation_height);
				}
				g.fillRect (0, y2, x2, line_separation_height);
			}
			else {
				g.fillRect (x1, y1, x2 - x1, line_separation_height);
			}
			if (y1 > h) break;
		}
		
		if (sequenceEditor != null) {
//			g.setColor (new Color (225, 174, 0));
			g.setColor (new Color (213, 165, 0));
//			g.setColor (new Color (0, 80, 0));
			int y1 = line_separation_offset + space * (sequenceEditor.getFirstIndex () / w) - offsetY * space;
			int x1 = sequenceEditor.getFirstIndex () % w;
			int y2 = line_separation_offset + space * (sequenceEditor.getLastIndex () / w) - offsetY * space;
			if (y2 >= 0 && y1 < h) {
				int x2 = sequenceEditor.getLastIndex () % w;
				if (y2 != y1) {
					g.fillRect (x1, y1 + line_separation_height / 2 - 2, w - x1, 4);
					for (int j = y1 + space; j < y2; j += space) {
						g.fillRect (0, j + line_separation_height / 2 - 2, w, 4);
					}
					g.fillRect (0, y2 + line_separation_height / 2 - 2, x2, 4);
				}
				else {
					g.fillRect (x1, y1 + line_separation_height / 2 - 2, x2 - x1, 4);
				}
				g.fillRect (x1, y1, 4, line_separation_height);
				g.fillRect (x2 - 4, y2, 4, line_separation_height);
			}
		}
		if (plotPairsJFrame.menuViewFeatureGene.isSelected ()) {
			for (int i = 0; i < reco.geneFeatures.size (); i++) {
				GeneFeature gene = reco.geneFeatures.get (i);
				int y1 = line_separation_offset + scaleY + space * (gene.start / w) - offsetY * space;
				int x1 = gene.start % w;
				int y2 = line_separation_offset + scaleY + space * (gene.stop / w) - offsetY * space;
				if (y2 < 0) continue;
				int x2 = gene.stop % w;
				
				String label = gene.gene;
				if (label != null) {
					int labelX;
					int labelY;
					if (y2 != y1) {
						if (y2 - y1 > space) {
							labelX = w / 2 - stringWidth (g, label) / 2;
							labelY = y1 + space;
						}
						else {
							if (w - x1 > x2) {
								labelY = y1;
								labelX = x1 + (w - x1) / 2 - stringWidth (g, label) / 2;
							}
							else {
								labelY = y2;
								labelX = x2 / 2 - stringWidth (g, label) / 2;
							}
						}
					}
					else {
						labelY = y1;
						labelX = x1 + (x2 - x1) / 2 - stringWidth (g, label) / 2;	
					}
					g.setColor (TRNA_AREA_COLOR);
					g.drawString (label, labelX, labelY - 2);
				}
				
				g.setColor (new Color (245, 245, 245));
				if (y2 != y1) {
					g.fillRect (x1, y1, w - x1, frames_height);
					for (int j = y1 + space; j < y2; j += space) {
						g.fillRect (0, j, w, frames_height);
					}
					g.fillRect (0, y2, x2, frames_height);
				}
				else {
					g.fillRect (x1, y1, x2 - x1, frames_height);
				}
				if (y1 > h) break;
				g.setColor (new Color (220, 220, 220));
				
				g.drawLine (x1, y1, x1, y1 + frames_height);
				g.drawLine (x2, y2, x2, y2 + frames_height);
				
			}
			g.setColor (new Color (200, 200, 200));
			for (int i = 0; i < reco.geneFeatures.size (); i++) {
				GeneFeature gene = reco.geneFeatures.get (i);
				int y1 = line_separation_offset + scaleY + space * (gene.start / w) - offsetY * space;
				int x1 = gene.start % w;
				int y2 = line_separation_offset + scaleY + space * (gene.stop / w) - offsetY * space;
				if (y2 < 0) continue;
				int x2 = gene.stop % w;
				int frame = gene.getFrame ();
				if (gene.reverse) {
					y2 += frame * frame_height + frame_height / 2;
					g.fillOval (x2 - 3, y2 - 3, 7, 7);
				} 
				else {
					y1 += frame * frame_height + frame_height / 2;
					g.fillOval (x1 - 3, y1 - 3, 7, 7);
				}
			}
		}
		
		g.setColor (FRAME_SEPARATION_COLOR);
		for (int i = 0; i < h; i += space) {
			for (int j = 1; j <= 6; j++) {
				g.drawLine (0, i + j * frame_height, w, i + j * frame_height);
			}
		}
		
		if (computeBox) orfBoxes.clear ();
		
		plotOrfs (g, reco.orfsPorFrame, w, h, space, frames_offset, frame_height);
		
		for (int fI = 0; fI < 6; fI ++) {
			List <Codon> codons = reco.startStopPorFrame [fI];
			if (plotPairsJFrame.menuViewStartStops.isSelected ()) {
				for (int i = 0; i < codons.size () - 1; i ++) {
					Codon c = codons.get (i);
					if (c.type == OrfRecognizer.VALINA || c.type == OrfRecognizer.VALINA_R) continue;
					int y1 = frames_offset + fI * frame_height + space * (c.index / w) - offsetY * space;
					if (y1 < 0) continue;
					int x1 = c.index % w;
					if (c.type == OrfRecognizer.START || c.type == OrfRecognizer.START_R)	g.setColor (COR_START_CODON);
	                else g.setColor (COR_STOP_CODON);
					g.drawLine (x1, y1, x1, y1 + scaleY);
					if (y1 > h) break;
				}
			}
			if (plotPairsJFrame.menuViewValina.isSelected ()) {
	            g.setColor (COR_VALINA);
	            for (int i = 0; i < codons.size () - 1; i ++) {
					Codon c = codons.get (i);
					if (c.type != OrfRecognizer.VALINA && c.type != OrfRecognizer.VALINA_R) continue;
					int y1 = frames_offset + fI * frame_height + space * (c.index / w) - offsetY * space;
					if (y1 < 0) continue;
					int x1 = c.index % w;
					g.drawLine (x1, y1, x1, y1 + scaleY);
					if (y1 > h) break;
				}
			}
		}
		
		if (getSelected () != null) {
            int ddx = 0;
            int ddy = 0;
            Font oldFont = g.getFont ();
			g.setStroke (SELECT_STROKE);
			g.setColor (COR_SELECT);
			g.setStroke (SELECT_STROKE);
			
			ORFBox bS = getSelected ();
			ORF orfS = bS.orf;
			
			
            for (ORFBox b: orfBoxes) {
                if (b.orf == orfS) g.fillRect (b.x, b.y, b.width + 1, b.height + 1);
            }
			g.setColor (COR_BCK_SELECT);
            String s = plotPairsJFrame.menuViewAcidoSequence.isSelected () ? orfS.getAminoSequence () : orfS.getSequence ();
            int size =  5 + (s.length () + 1) / 100;
            size = 5;
            
            if (bS.y - size * 12 - 4 < 0) ddy = bS.y + 4 + bS.height;
            else ddy = bS.y - size * 12 - 6;
            if (bS.x + 717 > w) ddx = w - 717;
            else ddx = bS.x;
            if (ddx < 2) ddx = 2;
                        
			g.fillRect (ddx, ddy, 715, size * 12 + 4);
            g.setColor (COR_SELECT);
            g.setStroke (NORMAL_STROKE);
            g.drawRect (ddx, ddy, 715, size * 12 + 4);
			g.setColor (Color.BLACK);
			g.drawString (orfS.id + " (length: " + s.length () + ")", ddx + 6, ddy + 12);
			g.drawString ("Entry: " + orfS.product.entryId + " Id: " + orfS.product.identity + " Sc.: " + orfS.product.score, ddx + 6, ddy + 24);
			g.drawString ("Name: " + orfS.product.productName, ddx + 6, ddy + 36);
			g.drawString ("Org: " + orfS.product.organism, ddx + 6, ddy + 60);
//			g.drawString ("Spe:" + orfS.getSpecificity (), ddx + 600, ddy +  48);
			int acc = (int) (orfS.getAccuracy () * 100);
			g.drawString ("Acc:" + acc / 100., ddx + 650, ddy +  36);
			if (orfS.product.hasMatchingSequence ()) {
				int boxW = 410;
				float [] boxes = orfS.getMatchBox ();
				g.fillRect ((int) (ddx + 300 + boxW * boxes [0]), ddy + 6, (int) (boxW * (boxes [1] - boxes [0])), 5);
				g.drawRect ((int) (ddx + 300 + boxW * boxes [4]), ddy + 3, (int) (boxW * (boxes [5] - boxes [4])), 17);
				float id = orfS.product.identity;
				if (id >= 50) {
					g.setColor (new Color ((int) (255 * id / 100), (int) (255 * (100 - id) / 50), 0));
				}
				else {
					g.setColor (new Color (0, (int) (255 * id / 50), (int) (255 * (100 - id) / 100)));
				}
				g.fillRect ((int) (ddx + 300 + boxW * boxes [2]), ddy + 13, (int) (boxW * (boxes [3] - boxes [2])), 5);
			}
			g.setColor (Color.BLACK);
			if (orfS.isGene ()) g.drawString ("Gene: " + orfS.product.gene, ddx + 6, ddy + 48);
			else g.drawString ("Product: " + orfS.product.productId, ddx + 6, ddy + 48);
            g.setFont (oldFont);
		}
		
		if (plotPairsJFrame.menuViewGC.isSelected ()) {
			g.setColor (new Color (0, 0, 128));
			int start = offsetY * w;
			int end = (offsetY + h / space) * w;
			if (mediaMovelGC == null) {
				mediaMovelGC = reco.getMediaGC (start, end);
			}
			int lastY = 0;
			int lastX = 0;
			for (int i = 0; i < mediaMovelGC.size (); i++) {
				float gc = mediaMovelGC.get (i) * 4 - 2;
				if (gc < 0) gc = 0;
				if (gc > 1) gc = 1;
				gc = (int) (frame_height * gc);
				int y = (int) gc + line_separation_offset + space * ((start + i) / w) - offsetY * space;
				int x = (start + i) % w;
				if (x > 0) g.drawLine (lastX, lastY, x, y);
				else g.drawLine (x, y, x, y);
				lastX = x;
				lastY = y;
			}
		}
		computeBox = false;
	}
	
	public void plotOrfs (Graphics2D g, List <ORF> orfsPFrame [], int w, int h, int space, int frames_offset, int frame_height) {
		boolean hideRemoved = plotPairsJFrame.menuViewHideRemovedOrf.isSelected ();
		
		for (int fI = 0; fI < 6; fI ++) {
			List <ORF> orfs = orfsPFrame [fI];
			for (int i = 0; i < orfs.size (); i ++) {
				ORF orf = orfs.get (i);
				if (orf.isRemoved () && hideRemoved) continue;
				
				int x1;
				int y1;
				int x2;
				int y2;
				
				orf.loadProduct ();
				y2 = frames_offset + fI * frame_height + space * (orf.stop / w) - offsetY * space;
				if (y2 < 0) continue;
				y1 = frames_offset + fI * frame_height + space * (orf.start / w) - offsetY * space;
				if (y1 > h) break;
					
				if (orf.product != null && orf.product.hasEntries) {
					float id = orf.product.identity;
					if (id >= 50) g.setColor (new Color ((int) (255 * id / 100), (int) (255 * (100 - id) / 50), 0));
					else g.setColor (new Color (0, (int) (255 * id / 50), (int) (255 * (100 - id) / 100)));
					
					int [] bound = orf.getQueryMatchBounds ();
					int startQuerySeq = bound [0];
					int endQuerySeq = bound [1];
					int startMatchSeq = bound [2];

					if (orf.isReverse ()) {
						y1 = frames_offset + fI * frame_height + space * ((orf.stop - endQuerySeq * 3) / w) - offsetY * space;
						x1 = (orf.stop - endQuerySeq * 3) % w;
						y2 = frames_offset + fI * frame_height + space * ((orf.stop - startQuerySeq * 3) / w) - offsetY * space;
						x2 = (orf.stop - startQuerySeq * 3) % w;
					}
					else {
						y1 = frames_offset + fI * frame_height + space * ((orf.start + startQuerySeq * 3) / w) - offsetY * space;
						x1 = (orf.start + startQuerySeq * 3) % w;
						y2 = frames_offset + fI * frame_height + space * ((orf.start + endQuerySeq * 3) / w) - offsetY * space;
						x2 = (orf.start + endQuerySeq * 3) % w;
					}
					if (y2 != y1) {
						g.fillRect (x1, y1 + 5, w - x1, frame_height - 10);
						for (int j = y1 + space; j < y2; j += space) g.fillRect (0, j + 5, w, frame_height - 10);
						g.fillRect (0, y2 + 5, x2, frame_height - 10);
					}
					else {
						g.fillRect (x1, y1 + 5, x2 - x1, frame_height - 10);
					}
					
					if (plotPairsJFrame.menuViewQueues.isSelected ()) {
						if (orf.isReverse()) {
							y1 = frames_offset + fI * frame_height + space * ((orf.stop - (startQuerySeq - startMatchSeq) * 3 - orf.product.matchingSequence.length () * 3) / w) - offsetY * space;
							x1 = (orf.stop - (startQuerySeq - startMatchSeq) * 3 - orf.product.matchingSequence.length () * 3) % w;
							y2 = frames_offset + fI * frame_height + space * ((orf.stop - (startQuerySeq - startMatchSeq) * 3) / w) - offsetY * space;
							x2 = (orf.stop - (startQuerySeq - startMatchSeq) * 3) % w;
						}
						else {
							y1 = frames_offset + fI * frame_height + space * ((orf.start + (startQuerySeq - startMatchSeq) * 3) / w) - offsetY * space;
							x1 = (orf.start + (startQuerySeq - startMatchSeq) * 3) % w;
							y2 = frames_offset + fI * frame_height + space * ((orf.start + (startQuerySeq - startMatchSeq) * 3 + orf.product.matchingSequence.length () * 3) / w) - offsetY * space;
							x2 = (orf.start + (startQuerySeq - startMatchSeq)* 3 + orf.product.matchingSequence.length () * 3) % w;
						}
						if (y2 != y1) {
							g.drawRect (x1, y1 + 5, w - x1, frame_height - 10);
							for (int j = y1 + space; j < y2; j += space) g.drawRect (0, j + 5, w, frame_height - 10);
							g.drawRect (0, y2 + 5, x2, frame_height - 10);
						}
						else {
							g.drawRect (x1, y1 + 5, x2 - x1, frame_height - 10);
						}
					}
				}
				
				if (fI >= 3) {
					if (orf.isRemoved ()) g.setColor (COR_REVERSE_WAY_REMOVED);
					else g.setColor (COR_REVERSE_WAY);
				}
				else {
					if (orf.isRemoved()) g.setColor (COR_NORMAL_WAY_REMOVED);
					else g.setColor (COR_NORMAL_WAY);
				}
				
				x1 = orf.start % w;
				y2 = frames_offset + fI * frame_height + space * (orf.stop / w) - offsetY * space;
				x2 = orf.stop % w;
				y1 = frames_offset + fI * frame_height + space * (orf.start / w) - offsetY * space;
				
				orf.loadProduct ();
				if (orf.product != null && orf.product.entryId != null) {
					boolean duplicated = reco.duplicatedEntries.isDuplicated (orf.product.entryId) || reco.duplicatedGenes.isDuplicated (orf.product.productId) ;
					Font ft = g.getFont ();
					if (orf.isGene ()) g.setFont (ft.deriveFont (Collections.singletonMap (TextAttribute.UNDERLINE, TextAttribute.UNDERLINE_ON)));
					else g.setFont (ft.deriveFont (Collections.singletonMap (TextAttribute.WEIGHT, TextAttribute.WEIGHT_EXTRA_LIGHT)));
					if (y2 != y1) {
						if (computeBox && (!orf.isRemoved () || !hideRemoved)) orfBoxes.add (new ORFBox (orf, x1, y1, w - x1, frame_height));
						if (!orf.product.blasted) g.fillRect (x1, y1, w - x1, frame_height);
						else g.drawRect (x1, y1, w - x1, frame_height);
						for (int j = y1 + space; j < y2; j += space) {	
							if (computeBox && (!orf.isRemoved () || !hideRemoved)) orfBoxes.add (new ORFBox (orf, 0, j, w, frame_height));
							if (!orf.product.blasted) g.fillRect (0, j, w, frame_height);
							else g.drawRect (0, j, w, frame_height);
						}
						if (computeBox && (!orf.isRemoved () || !hideRemoved)) orfBoxes.add (new ORFBox (orf, 0, y2, x2, frame_height));
						if (!orf.product.blasted) g.fillRect (0, y2, x2, frame_height);
						else g.drawRect (0, y2, x2, frame_height);
						if (!orf.isRemoved() && !orf.product.productId.equals ("")) {
							String label;
							if (orf.isGene ()) label = orf.product.gene + (duplicated ? " * " : "");
							else if (!orf.isUncharacterized ()) label = orf.product.productId + (duplicated ? " * " : "");
							else label = "?";
							if (orf.product.fragment) label = '~' + label;
							if (orf.product.swiss) label += " (S)";
							int labelX;
							int labelY;
							if (y2 - y1 > space) {
								labelX = w / 2 - stringWidth (g, label) / 2;
								labelY = y1 + space;
							}
							else {
								if (w - x1 > x2) {
									labelY = y1;
									labelX = x1 + (w - x1) / 2 - stringWidth (g, label) / 2;
								}
								else {
									labelY = y2;
									labelX = x2 / 2 - stringWidth (g, label) / 2;
								}
							}
							g.setColor (COR_SELECT);
							g.fillRect (labelX - 3, labelY + 1, stringWidth (g, label) + 6, 12);
							g.setColor (Color.BLACK);
							g.drawString (label, labelX, labelY + frame_height - 2);
						}
					}
					else {
						if (computeBox && (!orf.isRemoved () || !hideRemoved)) orfBoxes.add (new ORFBox (orf, x1, y1, x2 - x1, frame_height));
						if (!orf.product.blasted) g.fillRect (x1, y1, x2 - x1, frame_height);
						else g.drawRect (x1, y1, x2 - x1, frame_height);
						if (!orf.isRemoved() && !orf.product.productId.equals ("")) {
							String label;
							if (orf.isGene ()) label = orf.product.gene + (duplicated ? " * " : "");
							else if (!orf.isUncharacterized ()) label = orf.product.productId + (duplicated ? " * " : "");
							else label = "?";
							if (orf.product.fragment) label = '~' + label;
							if (orf.product.swiss) label += " (S)";
							int labelX = x1 + (x2 - x1) / 2 - stringWidth (g, label) / 2;
							g.setColor (COR_SELECT);
							g.fillRect (labelX - 3, y1 + 1, stringWidth (g, label) + 6, 12);
							g.setColor (Color.BLACK);
							g.drawString (label, labelX, y1 + frame_height - 2);
						}
					}
					g.setFont (ft);
				}
				if (fI >= 3) {
					Polygon pol = new Polygon ();
					pol.addPoint (x1, y1 + frame_height / 2);
					pol.addPoint (x1 + 10, y1);
					pol.addPoint (x1 + 10, y1 + frame_height);
					g.fillPolygon (pol);
					g.fillRect (x1, y1, 4, frame_height);
				}
				else {
					Polygon pol = new Polygon ();
					pol.addPoint (x2, y2 + frame_height / 2);
					pol.addPoint (x2 - 10, y2);
					pol.addPoint (x2 - 10, y2 + frame_height);
					g.fillPolygon (pol);
					g.fillRect (x2 - 3, y2, 3, frame_height);
				}
					
				if (orf.isRemoved () && orf.length >= CodonConfiguration.min_orf_length) {
					if (fI >= 3) g.setColor (COR_REVERSE_WAY);
					else g.setColor (COR_NORMAL_WAY);
					if (y2 != y1) {
						g.drawLine (x1, y1 + scaleY / 2, w, y1 + frame_height / 2);
						for (int j = y1 + space; j < y2; j += space) g.drawLine (0, j + frame_height / 2, w, j + frame_height / 2);
						g.drawLine (0, y2 + frame_height / 2, x2, y2 + frame_height / 2);
					}
					else {
						g.drawLine (x1, y1 + frame_height / 2, x2, y1 + frame_height / 2);
					}
				}
			}
		}
	}
	
	public static int stringWidth (Graphics2D g, String txt) {
		FontMetrics fMetric = g.getFontMetrics (g.getFont ());
		return fMetric.stringWidth (txt);
	}

	public void selectOrf (ORF orf) {
		int offsetY = orf.start / getWidth ();
		setOffsetY (offsetY, true);
		SwingUtilities.invokeLater (new Runnable () {
			public void run () {
				computeBox = true;
				repaint ();
				for (ORFBox box: orfBoxes) {
					if (box.orf == orf) {
						setSelected (box);
						break;
					}
				}
				repaint ();
			}
		});
	}

	public ORFBox getSelected () {
		return selected;
	}

	public void setSelected (ORFBox selected) {
		plotPairsJFrame.menuSelection.setEnabled (selected != null);
		this.selected = selected;
		if (selected != null) {
			plotPairsJFrame.menuSelectionForceRemoved.setSelected (selected.orf.isExcluded ());
			plotPairsJFrame.menuSelectionForceIncluded.setSelected (selected.orf.isForced ());
			sequenceEditor.setFirstIndex (selected.orf.start);
			sequenceEditor.repaint ();
		}
		
	}
	
}

