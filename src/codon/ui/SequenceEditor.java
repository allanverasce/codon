package codon.ui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.Toolkit;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.util.ArrayList;

import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;

import codon.data.ORF;
import codon.recognizer.OrfRecognizer;

public class SequenceEditor extends JPanel {
	private static final long serialVersionUID = 1L;
	
	OrfRecognizer reco;
	PlotOrfsFrame frame;
	
	private int firstIndex = 0;
	private int lastIndex = 0;
	
	EditorPainter painter;
	
	public SequenceEditor (PlotOrfsFrame frame) {
		this.frame = frame;
		this.reco = frame.reco;
		frame.plotOrfsPanel.sequenceEditor = this;
		setPreferredSize (new Dimension (500, 118));
		setLayout(new BorderLayout ());
		JButton bPrev = new JButton (new ImageIcon ("data/img/prevA.png"));
		JButton bNext = new JButton (new ImageIcon ("data/img/nextA.png"));
		painter = new EditorPainter ();
		add (bPrev, BorderLayout.WEST);
		add (painter, BorderLayout.CENTER);
		add (bNext, BorderLayout.EAST);
		bPrev.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				setFirstIndex (getFirstIndex () - 50);
				if (getFirstIndex () < 0) setFirstIndex (0);
				int lg = painter.getWidth () / 7 - 1;
				setLastIndex(getFirstIndex () + lg);
				frame.plotOrfsPanel.repaint ();
				painter.repaint ();
			}
		});
		bNext.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				setFirstIndex (getFirstIndex () + 50);
				int lg = painter.getWidth () / 7 - 1;
				if (getLastIndex () > reco.length ()) setFirstIndex (getLastIndex () - lg);
				frame.plotOrfsPanel.repaint ();
				painter.repaint ();
			}
		});
		addMouseWheelListener (new MouseWheelListener () {
			public void mouseWheelMoved (MouseWheelEvent e) {
				setFirstIndex (getFirstIndex () + e.getWheelRotation () * 20);
				if (getFirstIndex () < 0) setFirstIndex (0);
				int lg = painter.getWidth () / 7 - 1;
				if (getLastIndex () > reco.length ()) setFirstIndex (getLastIndex () - lg);
				frame.plotOrfsPanel.repaint ();
				painter.repaint ();
			}
		});
		setBorder (BorderFactory.createLineBorder (Color.BLACK));
		addComponentListener (new ComponentAdapter () {
			public void componentResized (ComponentEvent arg0) {
				int lg = painter.getWidth () / 7 - 1;
				setLastIndex (getFirstIndex () + lg);
				if (getLastIndex() > reco.length ()) setFirstIndex(getLastIndex() - lg);
				frame.plotOrfsPanel.repaint ();
			}
		});
	}
	
	public int getFirstIndex () {
		return firstIndex;
	}

	public void setFirstIndex (int firstIndex) {
		this.firstIndex = firstIndex;
		int lg = painter.getWidth () / 7 - 1;
		setLastIndex (getFirstIndex () + lg);
	}

	public int getLastIndex () {
		return lastIndex;
	}

	public void setLastIndex (int lastIndex) {
		this.lastIndex = lastIndex;
	}

	final static Color DARK_GREEN = new Color (0, 180, 0);
	final static Color DARK_BLUE = new Color (0, 0, 180);
	
	class EditorPainter extends JPanel {
		int selBegin = -1;
		int selEnd = -1;
		public EditorPainter () {
			setFocusable (true);
			JPopupMenu popup = new JPopupMenu ();
			JMenuItem menuEdit = new JMenuItem ("Edit");
			JMenuItem menuCopyForward = new JMenuItem ("Copy forward");
			JMenuItem menuCopyReverse = new JMenuItem ("Copy reverse");
			popup.add (menuEdit);
			popup.add (menuCopyForward);
			popup.add (menuCopyReverse);
			menuEdit.addActionListener (new ActionListener () {
				public void actionPerformed(ActionEvent e) {
					edit ();
				}
			});
			menuCopyForward.addActionListener (new ActionListener () {
				public void actionPerformed (ActionEvent e) {
					copyToClipBoard (reco.wholeSequence.substring (selBegin, selEnd));
				}
			});
			menuCopyReverse.addActionListener (new ActionListener () {
				public void actionPerformed (ActionEvent e) {
					copyToClipBoard (OrfRecognizer.reverseComplementar (reco.wholeSequence.substring (selBegin, selEnd)));
				}
			});
			addMouseListener (new MouseAdapter () {
				public void mouseClicked (MouseEvent e) {
					if (e.getButton () == MouseEvent.BUTTON3) {
						if (selBegin != selEnd) popup.show (EditorPainter.this, e.getX (), e.getY ());
					}
					if (e .getButton () == MouseEvent.BUTTON1 && e.getY () > 44 && e.getY () < 72) {
						if (e.getButton () == MouseEvent.BUTTON1 && !e.isShiftDown ()) {
							selBegin = firstIndex + (e.getX () - 2) / 7;
							selEnd = selBegin; 
						}
						if (e.getButton () == MouseEvent.BUTTON1 && e.isShiftDown ()) {
							selEnd = firstIndex + (e.getX () - 2) / 7 + 1;
//							
						}
						repaint ();
						grabFocus ();
					}
				}
			});
			addKeyListener (new KeyAdapter () {
				public void keyPressed (KeyEvent e) {
					if (e.isControlDown () && e.getKeyCode () == KeyEvent.VK_E) {
						edit ();
					}
					if (e.isControlDown () && e.getKeyCode () == KeyEvent.VK_C) {
						copyToClipBoard (reco.wholeSequence.substring (selBegin, selEnd));
					}
					if (e.isControlDown () && e.getKeyCode () == KeyEvent.VK_R) {
						copyToClipBoard (OrfRecognizer.reverseComplementar (reco.wholeSequence.substring (selBegin, selEnd)));
					}
				}
				public void keyReleased (KeyEvent e) {
					if (e.getKeyCode () == KeyEvent.VK_SHIFT) {
						if (selBegin > selEnd) {
							int tmp = selEnd;
							selEnd = selBegin;
							selBegin = tmp;
						}
					}
				}
			});
			
		}
		
		protected void copyToClipBoard (String s) {
			 Toolkit.getDefaultToolkit ().getSystemClipboard ().setContents (new StringSelection (s), null);;
		}

		public void edit () {
			if (selBegin != selEnd) {
				String newS = JOptionPane.showInputDialog (PlotOrfsFrame.instance, "Input the new sequence that will substitute the actual", reco.wholeSequence.substring (selBegin, selEnd));
				if (newS == null) return;
				newS = newS.toUpperCase ();
				for (int i = 0; i < newS.length (); i++) {
					switch (newS.charAt (i)) {
					case 'A':
					case 'T':
					case 'C':
					case 'G':
						break;
					default:
						JOptionPane.showMessageDialog (PlotOrfsFrame.instance, "Illegal character " + newS.charAt (i) + " the alteration will be ignored	");
						return;
					}
				}
				if (newS != null) {
					reco.replace (selBegin, selEnd, newS);
					selEnd = selBegin + newS.length ();
					repaint ();
					frame.plotOrfsPanel.repaint ();
				}
			}
		}
		
		private static final long serialVersionUID = 1L;
		Font fontMono = new Font (Font.MONOSPACED, Font.PLAIN, 12);
		protected void paintComponent (Graphics g) {
			super.paintComponent (g);
			((Graphics2D) g).setRenderingHint (RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
			g.setFont (fontMono);
			int w = getWidth ();
			int lg = w / 7 - 1;
			lastIndex = Math.min (lastIndex, reco.length ());
			
			g.setColor (Color.WHITE);
			g.fillRect (0, 44, w, 14);
			g.setColor (Color.DARK_GRAY);
			g.fillRect (0, 58, w, 14);
			
			g.setColor (Color.LIGHT_GRAY);
			int sB = selBegin;
			int sE = selEnd;
			if (sB > sE) {
				int tmp = sE;
				sE = sB;
				sB = tmp;
			}
			int indexStart = (sB - getFirstIndex ());
			int posStart = indexStart * 7;
			int indexStop = (sE - getFirstIndex ());
			int posStop =  indexStop * 7;
			g.setColor (Color.LIGHT_GRAY);
			g.fillRect (posStart + 2, 44, posStop - posStart, 28);
			g.setColor (Color.DARK_GRAY);
			g.drawRect (posStart + 2, 44, posStop - posStart, 28);
			
			int lgMax = reco.wholeSequence.length ();
			g.setColor (Color.BLACK);
			g.drawString (reco.wholeSequence.substring (Math.min (lgMax, getFirstIndex ()), Math.min (lgMax, lastIndex)), 2, 56);
			g.setColor (Color.WHITE);
			g.drawString (OrfRecognizer.complementar (reco.wholeSequence.substring (Math.min (lgMax, getFirstIndex ()), Math.min (lgMax, lastIndex))), 2, 70);
			g.setColor (Color.BLACK);
			
			for (ArrayList <ORF> orfs: reco.orfsPorFrame) {
				for (ORF orf: orfs) {
					if (orf.isRemoved ()) continue;
					if (orf.start >= getFirstIndex () && orf.start <= lastIndex || orf.stop >= getFirstIndex () && orf.stop <= lastIndex || orf.start <= getFirstIndex () && orf.stop >= lastIndex) {
						float id = 0;
						if (orf.product != null) id = orf.product.identity;
						if (id >= 50) g.setColor (new Color ((int) (255 * id / 100), (int) (255 * (100 - id) / 50), 0));
						else g.setColor (new Color (0, (int) (255 * id / 50), (int) (255 * (100 - id) / 100)));
						indexStart = (orf.start - getFirstIndex ());
						posStart = indexStart * 7;
						indexStop = (orf.stop - getFirstIndex ()) + 3; 
						posStop =  indexStop * 7;
						if (orf.frame < 3) g.fillRect (posStart, 3 + orf.frame * 14, posStop - posStart, 13);
						if (orf.frame >= 3) g.fillRect (posStart, 31 + orf.frame * 14, posStop - posStart, 13);
					}
					if (orf.start > lastIndex) break;
				}
			}
			
			for (int i = 0; i < lg - 1; i += 3) {
				String bases0 = reco.wholeSequence.substring (Math.min (lgMax, getFirstIndex () + i), Math.min (lgMax,getFirstIndex () + i + 3));
				String bases1 = reco.wholeSequence.substring (Math.min (lgMax, getFirstIndex () + i + 1), Math.min (lgMax, getFirstIndex () + i + 4));
				String bases2 = reco.wholeSequence.substring (Math.min (lgMax, getFirstIndex () + i + 2), Math.min (lgMax, getFirstIndex () + i + 5));
				
				int f = getFirstIndex ();
				if (getFirstIndex () < 2) f += 3;
				int frame0 = (f - 2) % 3;
				int frame1 = (f + 1 - 2) % 3;
				int frame2 = (f + 2 - 2) % 3;
				
				String amino0 = OrfRecognizer.baseAcidoMap.get (bases0);
				String amino1 = OrfRecognizer.baseAcidoMap.get (bases1);
				String amino2 = OrfRecognizer.baseAcidoMap.get (bases2);
				String amino3 = OrfRecognizer.baseAcidoMap.get (OrfRecognizer.reverseComplementar (bases0));
				String amino4 = OrfRecognizer.baseAcidoMap.get (OrfRecognizer.reverseComplementar (bases1));
				String amino5 = OrfRecognizer.baseAcidoMap.get (OrfRecognizer.reverseComplementar (bases2));
				
				g.setColor (getColor (amino0));
				if (amino0 != null) g.drawString (amino0, i * 7 + 7 + 2, 14 + frame0 * 14);
				g.setColor (getColor (amino1));
				if (amino1 != null) g.drawString (amino1, i * 7 + 14 + 2, 14 + frame1 * 14);
				g.setColor (getColor (amino2));
				if (amino2 != null) g.drawString (amino2, i * 7 + 21 + 2, 14 + frame2 * 14);
				g.setColor (getColor (amino3));
				if (amino3 != null) g.drawString (amino3, i * 7 + 7 + 2, 84 + frame0 * 14);
				g.setColor (getColor (amino4));
				if (amino4 != null) g.drawString (amino4, i * 7 + 14 + 2, 84 + frame1 * 14);
				g.setColor (getColor (amino5));
				if (amino5	 != null) g.drawString (amino5, i * 7 + 21 + 2, 84 + frame2 * 14);
				g.setColor (Color.BLACK);
			}
		}
		
		public Color getColor (String am) {
			if (am == null) return Color.BLACK;	
			switch (am.charAt (0)) {
			case '+':
			case '*':
			case '#':
				return DARK_BLUE;
			case 'M':
				return DARK_GREEN;
			case 'V':
			case 'L':
				return Color.YELLOW;
			default:
				return Color.BLACK;
			}
		}
	}
	
	
}

