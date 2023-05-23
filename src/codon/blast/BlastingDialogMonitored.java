package codon.blast;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import javax.sound.sampled.AudioFormat;
import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;
import javax.sound.sampled.Clip;
import javax.sound.sampled.DataLine;
import javax.sound.sampled.LineEvent;
import javax.sound.sampled.LineListener;
import javax.sound.sampled.LineUnavailableException;
import javax.sound.sampled.UnsupportedAudioFileException;
import javax.swing.BorderFactory;
import javax.swing.DefaultListCellRenderer;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.ScrollPaneConstants;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.Timer;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import codon.data.CodonConfiguration;
import codon.data.ORF;
import codon.recognizer.OrfRecognizer;
import codon.ui.PlotOrfsFrame;

public class BlastingDialogMonitored extends JPanel implements BlastingMonitorListener {
	private static final long serialVersionUID = 1L;
	
	BlastingMonitor blastingMonitor = null;
	
	class DoneBlast {
		ORF orf;
		String time;
		boolean visited = false;
		public DoneBlast (ORF orf, long time) {
			this.orf = orf;
			String sec = time + "";
			while (sec.length () < 4) sec = " " + sec;
			this.time = sec;
		}
		public String toString () {
			return time + "s - " + orf.id;
		}
	}
	
	class ProcessingBlast {
		ORF orf;
		long startTime;
		public boolean status;
		int tries = 1;
		public ProcessingBlast (ORF orf) {
			this.orf = orf;
			startTime = System.currentTimeMillis ();
		}
		public String toString () {
			long tm = System.currentTimeMillis () - startTime;
			tm /= 1000;
			return tm + "s " + (tries == 1 ? "": "(" + tries + " tries) ") + "- " + (status ? "(paused) - " : "") + orf.id;	
		}
	}
	
	Vector <ProcessingBlast> processingBlasts = new Vector <> ();
	
	class DoneBlastListRenderer extends DefaultListCellRenderer {
		private static final long serialVersionUID = 1L;
		public Component getListCellRendererComponent (JList <?> list, Object value, int index, boolean isSelected, boolean cellHasFocus) {
			Component cmp = super.getListCellRendererComponent (list, value, index, isSelected, cellHasFocus);
			JPanel p = new JPanel (new BorderLayout ());
			p.add (cmp);
			DoneBlast dB = (DoneBlast) value;
			JLabel l = new JLabel (dB.orf.isUncharacterized () ? " ? " : (dB.orf.isGene () ? "   " : " P "));
			l.setOpaque(true);
			p.add (l, BorderLayout.WEST);
			float accuracy = dB.orf.getAccuracy ();
			accuracy = Math.max (0, accuracy);
			accuracy = Math.min (100, accuracy);
			if (accuracy >= 50) l.setBackground (new Color ((int) (255 * accuracy / 100), (int) (255 * (100 - accuracy) / 50), 0));
			else l.setBackground (new Color (0, (int) (255 * accuracy / 50), (int) (255 * (100 - accuracy) / 100)));
			if (!dB.orf.isBlasted () || !dB.orf.product.hasEntries) {
				cmp.setBackground (Color.RED);
			}
			else if (!dB.visited) cmp.setBackground (Color.LIGHT_GRAY);
			else cmp.setBackground (Color.WHITE);
			return p;
		}
	}
	
	Vector <DoneBlast> doneOrfs = new Vector <> ();
	
	JLabel lToProcessNb = new JLabel ();
	JLabel lIgnoredNb = new JLabel ();
	JLabel lDoneNb = new JLabel ();
	JLabel lPreviouslyDoneNb = new JLabel ();
	PlotOrfsFrame plotOrfsFrame;
	JButton buttonCancel = new JButton ("Cancel");
	JList <ProcessingBlast> processingList = new JList <> ();
	JList <DoneBlast> doneList = new JList <> ();
	
	boolean changed = false;
	
	OrfRecognizer reco;
	JSpinner spiThread = new JSpinner ();
	
	JCheckBox ckPlaySound = new JCheckBox ("Play sound when blast ended");
	JCheckBox ckSpeedUpWithSwiss = new JCheckBox ("Speed-up with SwissProt");
	
	BlastHistoryPanel historyPanel = new BlastHistoryPanel ();
	boolean fullBlast;
	
	public BlastingDialogMonitored (PlotOrfsFrame plotOrfsFrame, List <ORF> validOrfs, int threadNb, boolean accelerateSwissProt, boolean fullBlast) {
		this.reco = plotOrfsFrame.reco;
		this.plotOrfsFrame = plotOrfsFrame;
		
		blastingMonitor = new BlastingMonitor (reco.dataSrcFileS, validOrfs.size (), threadNb, accelerateSwissProt, fullBlast, true, -1, -1);
		blastingMonitor.addBlastingMonitorListener (this);
		
		setLayout (new BorderLayout ());
		
		JPanel statsPanel = new JPanel (new GridLayout (6, 2));
		statsPanel.add (new JLabel (" Simultaneous tries"));
		statsPanel.add (spiThread);
		spiThread.setPreferredSize (new Dimension (100, 22));
		statsPanel.add (new JLabel (" To process"));
		statsPanel.add (lToProcessNb);
		statsPanel.add (new JLabel (" Discarded"));
		statsPanel.add (lIgnoredNb);
		statsPanel.add (new JLabel (" Done"));
		statsPanel.add (lDoneNb);
		statsPanel.add (new JLabel (" Previously Done"));
		statsPanel.add (lPreviouslyDoneNb);
		
		SpinnerModel sm = new SpinnerNumberModel (threadNb, 1, 20, 1);
		spiThread.setModel (sm);
		
		JPanel listPanel = new JPanel (new GridLayout (2, 1));
		JPanel pCurrent = new JPanel (new BorderLayout ());
		pCurrent.setBorder (BorderFactory.createTitledBorder ("Current blasting"));
//		pCurrent.add (new JLabel ("Current blasting"));
		JScrollPane scroll = new JScrollPane (processingList);
		pCurrent.add (scroll);
		JPanel pEnded = new JPanel (new BorderLayout ());
//		pEnded.add (new JLabel ("Ended"));
		pEnded.setBorder (BorderFactory.createTitledBorder ("Ended"));
		scroll = new JScrollPane (doneList);
		doneList.setCellRenderer (new DoneBlastListRenderer ());
		pEnded.add (scroll);
		JPanel panelCk = new JPanel (new GridLayout (2, 1));
		panelCk.add (ckPlaySound);
		if (CodonConfiguration.ENABLE_SWISS) panelCk.add (ckSpeedUpWithSwiss);
		ckSpeedUpWithSwiss.setSelected (accelerateSwissProt);
		pEnded.add (panelCk, BorderLayout.NORTH);
		ckPlaySound.setSelected (false);
		
		scroll = new JScrollPane (historyPanel);
		scroll.setHorizontalScrollBarPolicy (ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
		scroll.setVerticalScrollBarPolicy (ScrollPaneConstants.VERTICAL_SCROLLBAR_NEVER);
		pEnded.add (scroll, BorderLayout.SOUTH);
		
		listPanel.add (pCurrent);
		listPanel.add (pEnded);
		
		add (statsPanel, BorderLayout.NORTH);
		add (listPanel);
		
		JPopupMenu popupThreadList = new JPopupMenu ();
		JMenuItem menuKillItem = new JMenuItem ("Stop");
		JMenuItem menuKillAll = new JMenuItem ("Restart All");
		popupThreadList.add (menuKillItem);
		popupThreadList.add (menuKillAll);
		
		processingList.addMouseListener (new MouseAdapter () {
			public void mouseClicked (MouseEvent e) {
				if (e.getButton () == MouseEvent.BUTTON1 && e.getClickCount () > 1 && processingList.getSelectedIndex () != -1) {
					plotOrfsFrame.plotOrfsPanel.selectOrf (processingList.getSelectedValue ().orf);
				}
				if (e.getButton () == MouseEvent.BUTTON3) {
					popupThreadList.show (processingList, e.getX (), e.getY ());
				}
			}
		});
		processingList.addKeyListener (new KeyAdapter () {
			public void keyPressed (KeyEvent e) {
				if (e.getKeyCode () == KeyEvent.VK_DELETE) {
					ProcessingBlast p = processingList.getSelectedValue ();
					if (p != null) {
						blastingMonitor.restartOrf (p.orf.id);
						processingBlasts.remove (processingList.getSelectedIndex ());
						processingList.setListData (processingBlasts);
					}
				}
			}
		});
		menuKillItem.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				ProcessingBlast p = processingList.getSelectedValue ();
				if (p != null) {
					blastingMonitor.restartOrf (p.orf.id);
					processingBlasts.remove (processingList.getSelectedIndex ());
					processingList.setListData (processingBlasts);
				}
			}
		});
		menuKillAll.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				blastingMonitor.restartAll ();
			}
		});
		doneList.addMouseListener (new MouseAdapter () {
			public void mouseClicked (MouseEvent e) {
				if (e.getClickCount () > 1 && doneList.getSelectedIndex () != -1) {
					doneOrfs.get (doneList.getSelectedIndex ()).visited = true;
					plotOrfsFrame.plotOrfsPanel.selectOrf (doneOrfs.get (doneList.getSelectedIndex ()).orf);
					doneList.repaint ();
				}
			}
		});
		doneList.addKeyListener(new KeyAdapter () {
			public void keyPressed (KeyEvent e) {
				if (e.getKeyCode () == KeyEvent.VK_SPACE) {
					doneOrfs.get (doneList.getSelectedIndex ()).visited = true;
				}
				else if (e.getKeyCode () == KeyEvent.VK_ENTER) {
					doneOrfs.get (doneList.getSelectedIndex ()).visited = true;
					plotOrfsFrame.plotOrfsPanel.selectOrf (doneOrfs.get (doneList.getSelectedIndex ()).orf);
					doneList.repaint ();
				}
			}
		});
		ckSpeedUpWithSwiss.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				blastingMonitor.setAccelerateSwiss (ckSpeedUpWithSwiss.isSelected ());
			}
		});
		sm.addChangeListener (new ChangeListener () {
			public void stateChanged (ChangeEvent e) {
				int nb = (int) sm.getValue ();
				blastingMonitor.setSimultaneousTries (nb);
			}
		});
		Timer tm = new Timer (1000, new ActionListener () {
			public void actionPerformed (ActionEvent e) {
				ckSpeedUpWithSwiss.setSelected (blastingMonitor.isAccelerateSwiss ());
				if (changed) {
					reco.computeOrfFreeAndOrfSuperpositionSpans (plotOrfsFrame.tabRight.overlapPane.getTolerance ()); 
			    	processingList.setListData (processingBlasts);
			    	doneList.setListData (doneOrfs);
			    	playEndSound ();
			    	plotOrfsFrame.plotOrfsPanel.repaint ();
					updateStats ();
					changed = false;
				}
				else processingList.updateUI ();
			}
		});
		tm.setRepeats (true);
		tm.start ();
		add (buttonCancel, BorderLayout.SOUTH);
		buttonCancel.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				blastingMonitor.cancel ();
				plotOrfsFrame.tabRight.removeBlastTab (BlastingDialogMonitored.this);
				setVisible (false);
			}
		});
		plotOrfsFrame.tabRight.addBlastTab (this);
		blastingMonitor.start ();
	}
	
	public void playEndSound () {
		if (ckPlaySound.isSelected ()) {
			new Thread () {
				boolean playCompleted = false; 
				public void run () {
					try {
		  				File audioFile = new File ("data/end.wav");
		  				AudioInputStream audioStream;
		  				audioStream = AudioSystem.getAudioInputStream (audioFile);
		  				AudioFormat format = audioStream.getFormat ();
		  				DataLine.Info info = new DataLine.Info (Clip.class, format);
		  				Clip audioClip = (Clip) AudioSystem.getLine (info);
		  				audioClip.open (audioStream);
		  				audioClip.addLineListener (new LineListener () {
		  					public void update (LineEvent event) {
		  				        if (event.getType () == LineEvent.Type.STOP) {
		  				            playCompleted = true;
		  				        }
		  				    }
		  				});
  			  			audioClip.start ();
  		  				 while (!playCompleted) {
  		  	                try {
  		  	                    Thread.sleep (50);
  		  	                } catch (InterruptedException ex) {
  		  	                    ex.printStackTrace ();
  		  	                }
  		  	            }
  		  				audioClip.close ();
  		  				audioStream.close ();
  		  			} catch (UnsupportedAudioFileException e) {
  		  				e.printStackTrace ();
  		  			} catch (IOException e) {
  		  				e.printStackTrace ();
  		  			} catch (LineUnavailableException e) {
  		  			}
				}
			}.start ();
  		}
	}
		
	public void updateStats () {
		lToProcessNb.setText (blastingMonitor.toProcess + "");
		lIgnoredNb.setText (blastingMonitor.ignored + "");
		lDoneNb.setText (blastingMonitor.ended + "");
		lPreviouslyDoneNb.setText (blastingMonitor.previouslyDone + "");
	}
	
	public void blastStarted (BlastingMonitor manager, String orfId) {
		ORF blasted = null;
		for (ArrayList <ORF> orfs: reco.orfsPorFrame) for (ORF orf: orfs) {
			if (orf.id.equals (orfId)) {
				blasted = orf;
				break;
			}
		}
		if (blasted == null) {
			System.err.println ("Unknown ORF " + orfId);
		}
		else {
			processingBlasts.add (new ProcessingBlast (blasted));
		}
//		SwingUtilities.invokeLater (new Runnable () {
//			public void run () {
//				reco.computeOrfFreeAndOrfSuperpositionSpans (plotOrfsFrame.tabRight.overlapPane.getTolerance ()); 
//		    	processingList.setListData (processingBlasts);
//		    	doneList.setListData (doneOrfs);
//		    	playEndSound ();
//		    	plotOrfsFrame.plotOrfsPanel.repaint ();
//				updateStats ();
//			}
//		});
		changed = true;
	}

	public void blastEnded (BlastingMonitor manager, String orfId) {
		historyPanel.newBlast ();
		ORF mainOrf = null;
		long startTime = 0;
		for (ProcessingBlast pB : processingBlasts) {
			if (pB.orf.id.equals (orfId)) {
				mainOrf = pB.orf;
				startTime = pB.startTime;
				processingBlasts.remove (pB);
				break;
			}
		}
		if (mainOrf == null) {
			for (ArrayList <ORF> orfs: reco.orfsPorFrame) for (ORF orf: orfs) {
				if (orf.id.equals (orfId)) {
					mainOrf = orf;
					break;
				}
			}
		}
		if (mainOrf == null) return;
		mainOrf.reloadProduct ();
		mainOrf.autoRebound ();
		doneOrfs.add (0, new DoneBlast (mainOrf, (System.currentTimeMillis () - startTime) / 1000));
		changed = true; 
	}

	public void blastingFinalized (BlastingMonitor manager) {
		updateStats ();
		plotOrfsFrame.updateFilter ();
		JOptionPane.showMessageDialog (plotOrfsFrame, "Blasting finalized");
	}

	public void orfIgnored (BlastingMonitor manager, String orfId) {
		for (ArrayList <ORF> orfs: reco.orfsPorFrame) for (ORF orf: orfs) {
			if (orf.id.equals (orfId)) {
				orf.setRemoved (true);
				changed = true;
				break;
			}
		}
	}

	public void blastPreviouslyDone (BlastingMonitor manager, String orfId) {
		for (ArrayList <ORF> orfs: reco.orfsPorFrame) for (ORF orf: orfs) {
			if (orf.id.equals (orfId)) {
				orf.reloadProduct ();
//				orf.selectBetterEntryMaximizingIntergenicRegion ();
				orf.autoRebound ();
				changed = true;
//				SwingUtilities.invokeLater (new Runnable () {
//					public void run () {
//						doneOrfs.add (0, new DoneBlast (orf, 0));
//						doneList.setListData (doneOrfs);
//						updateStats ();
//					}
//				});
				break;
			}
		}
	}
	
	public void orfRebounded (BlastingMonitor manager, String orfId, int start, int stop) {
		for (ArrayList <ORF> orfs: reco.orfsPorFrame) for (ORF orf: orfs) {
			if (orf.id.equals (orfId)) {
				orf.start = start;
				orf.stop = stop;
			}
		}
	}

	public void orfReboundedAndRenamed (BlastingMonitor manager, String orfId, int start, int stop) {
		for (ArrayList <ORF> orfs: reco.orfsPorFrame) for (ORF orf: orfs) {
			if (orf.id.equals (orfId)) {
				orf.setOriginalStart (start);
				orf.setOriginalStop (stop);
				orf.rename ();
			}
		}
	}

	public void blastRestarted (BlastingMonitor manager) {
//		doneOrfs.clear ();
		processingBlasts.clear ();
		changed = true;
//		SwingUtilities.invokeLater (new Runnable () {
//			public void run () {
//				
//				processingList.setListData (processingBlasts);
//		    	doneList.setListData (doneOrfs);
//		    	updateStats ();
//			}
//		});
	}

	public void orfReincluded (BlastingMonitor manager, String orfId) {
		for (ArrayList <ORF> orfs: reco.orfsPorFrame) for (ORF orf: orfs) {
			if (orf.id.equals (orfId)) {
				orf.setRemoved (false);
				changed = true;	
				break;
			}
		}
	}

	public void blastRestarted (BlastingMonitor manager, String orfId) {
		for (ProcessingBlast pB : processingBlasts) {
			if (pB.orf.id.equals (orfId)) {
				pB.startTime = System.currentTimeMillis ();
				pB.tries ++;
				break;
			}
		}
	}

	public void statusChanged (BlastingMonitor blastingMonitor, String orfId, boolean status) {
		for (ProcessingBlast pB : processingBlasts) {
			if (pB.orf.id.equals (orfId)) {
				pB.status = status;
				break;
			}
		}
	}
	
}
