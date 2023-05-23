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
import javax.swing.SwingUtilities;
import javax.swing.Timer;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import codon.data.ORF;
import codon.recognizer.OrfRecognizer;
import codon.ui.PlotOrfsFrame;

public class BlastingDialog extends JPanel implements BlastingManagerListener {
	private static final long serialVersionUID = 1L;
	
	BlastingManager blastingManager = null;
	
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
	JList <BlastingThread> threadsList = new JList <> ();
	JList <DoneBlast> doneList = new JList <> ();
	
	OrfRecognizer reco;
	JSpinner spiThread = new JSpinner ();
	
	JCheckBox ckPlaySound = new JCheckBox ("Play sound when blast ended");
	JCheckBox ckSpeedUpWithSwiss = new JCheckBox ("Speed-up with SwissProt");
	
	BlastHistoryPanel historyPanel = new BlastHistoryPanel ();
	boolean fullBlast;
	
	public BlastingDialog (PlotOrfsFrame plotOrfsFrame, List <ORF> validOrfs, int threadNb, boolean upgradeWithSwissProt, boolean fullBlast) {
		this (plotOrfsFrame, threadNb, upgradeWithSwissProt);
		this.fullBlast = fullBlast;
		blastingManager.fullBlast = fullBlast;
		add (buttonCancel, BorderLayout.SOUTH);
		buttonCancel.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				blastingManager.interrupt ();
				plotOrfsFrame.tabRight.removeBlastTab (BlastingDialog.this);
				setVisible (false);
			}
		});
		plotOrfsFrame.tabRight.addBlastTab (this);
		start (validOrfs);
	}

	public BlastingDialog (PlotOrfsFrame plotOrfsFrame, int threadNb, boolean accelerateSwissProt) {
		this.reco = plotOrfsFrame.reco;
		this.plotOrfsFrame = plotOrfsFrame;
		
		blastingManager = new BlastingManager (reco, threadNb, fullBlast);
		blastingManager.addBlastingManagerListener (this);
		
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
		
		SpinnerModel sm = new SpinnerNumberModel (threadNb, 1, 100, 1);
		spiThread.setModel (sm);
		
		JPanel listPanel = new JPanel (new GridLayout (2, 1));
		JPanel pCurrent = new JPanel (new BorderLayout ());
		pCurrent.setBorder (BorderFactory.createTitledBorder ("Current blasting"));
//		pCurrent.add (new JLabel ("Current blasting"));
		JScrollPane scroll = new JScrollPane (threadsList);
		pCurrent.add (scroll);
		JPanel pEnded = new JPanel (new BorderLayout ());
//		pEnded.add (new JLabel ("Ended"));
		pEnded.setBorder (BorderFactory.createTitledBorder ("Ended"));
		scroll = new JScrollPane (doneList);
		doneList.setCellRenderer (new DoneBlastListRenderer ());
		pEnded.add (scroll);
		JPanel panelCk = new JPanel (new GridLayout (2, 1));
		panelCk.add (ckPlaySound);
		panelCk.add (ckSpeedUpWithSwiss);
		pEnded.add (panelCk, BorderLayout.NORTH);
		ckPlaySound.setSelected (false);
		ckSpeedUpWithSwiss.setSelected (accelerateSwissProt);
		Blaster.acceleratedWithSwissProt = accelerateSwissProt;
		
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
		popupThreadList.add (menuKillItem);
		
		threadsList.addMouseListener (new MouseAdapter () {
			public void mouseClicked (MouseEvent e) {
				if (e.getButton () == MouseEvent.BUTTON1 && e.getClickCount () > 1 && threadsList.getSelectedIndex () != -1) {
					plotOrfsFrame.plotOrfsPanel.selectOrf (blastingManager.threads.get (threadsList.getSelectedIndex ()).orf);
				}
				if (e.getButton () == MouseEvent.BUTTON3) {
					popupThreadList.show (threadsList, e.getX (), e.getY ());
				}
			}
		});
		threadsList.addKeyListener (new KeyAdapter () {
			public void keyPressed (KeyEvent e) {
				if (e.getKeyCode () == KeyEvent.VK_DELETE) {
					int ind = threadsList.getSelectedIndex ();
					if (ind != -1) {
						blastingManager.kill (ind);
					}
				}
			}
		});
		menuKillItem.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				int ind = threadsList.getSelectedIndex ();
				if (ind != -1) {
					blastingManager.kill (ind);
				}
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
				Blaster.acceleratedWithSwissProt = ckSpeedUpWithSwiss.isSelected ();
			}
		});
		sm.addChangeListener (new ChangeListener () {
			public void stateChanged (ChangeEvent e) {
				int nb = (int) sm.getValue ();
				blastingManager.setMaxThread (nb);
			}
		});
		Timer tm = new Timer (1000, new ActionListener () {
			public void actionPerformed (ActionEvent e) {
				ckSpeedUpWithSwiss.setSelected (Blaster.acceleratedWithSwissProt);
				threadsList.updateUI ();
			}
		});
		tm.setRepeats (true);
		tm.start ();
	}
	
	public void start (List <ORF> validOrfs) {
		blastingManager.start (validOrfs);
		plotOrfsFrame.updateFilter ();
		updateStats ();
		threadsList.setListData (blastingManager.threads);
	}
	
	synchronized public void ended (BlastingThread thread) {
		historyPanel.newBlast ();
		ORF mainOrf = thread.orf; 
		doneOrfs.add (0, new DoneBlast (mainOrf, (System.currentTimeMillis () - thread.startTime) / 1000));
		plotOrfsFrame.plotOrfsPanel.repaint ();
		updateStats ();
		SwingUtilities.invokeLater (new Runnable () {
			public void run () {
				reco.computeOrfFreeAndOrfSuperpositionSpans (plotOrfsFrame.tabRight.overlapPane.getTolerance ()); 
		    	threadsList.setListData (blastingManager.threads);
		    	doneList.setListData (doneOrfs);
		    	threadsList.updateUI ();
		    	playEndSound ();
			}
		});
	}
	
	synchronized public void addOrfToBlast (ORF orf) {
		blastingManager.addOrfToBlast (orf);
		updateStats ();
		SwingUtilities.invokeLater (new Runnable () {
	      public void run () {
	    	  threadsList.setListData (blastingManager.threads);
	    	  threadsList.updateUI ();
	      }
	     });
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
		lToProcessNb.setText (blastingManager.toProcess () + "");
		lIgnoredNb.setText (blastingManager.ignored + "");
		lDoneNb.setText (blastingManager.done + "");
		lPreviouslyDoneNb.setText (blastingManager.previouslyDone + "");
	}

	public void threadStatustChanged (BlastingManager manager, BlastingThread endedThread) {
		SwingUtilities.invokeLater (new Runnable () {
			public void run () {
				threadsList.repaint ();						
			}
		});
	}

	public void threadEnded (BlastingManager manager, BlastingThread endedThread) {
		historyPanel.newBlast ();
		ORF mainOrf = endedThread.orf; 
		doneOrfs.add (0, new DoneBlast (mainOrf, (System.currentTimeMillis () - endedThread.startTime) / 1000));
		plotOrfsFrame.plotOrfsPanel.repaint ();
		updateStats ();
		SwingUtilities.invokeLater (new Runnable () {
			public void run () {
				reco.computeOrfFreeAndOrfSuperpositionSpans (plotOrfsFrame.tabRight.overlapPane.getTolerance ()); 
		    	threadsList.setListData (blastingManager.threads);
		    	doneList.setListData (doneOrfs);
		    	threadsList.updateUI ();
		    	playEndSound ();
			}
		});
	}

	public void blastingFinalized (BlastingManager manager) {
		JOptionPane.showMessageDialog (plotOrfsFrame, "Blasting finalized");
	}

	public void threadStarted (BlastingManager manager, BlastingThread endedThread) {}
	public void threadRestarted (BlastingManager manager, BlastingThread endedThread) {}
	public void orfIgnored (BlastingManager manager, ORF orf) {}
	public void orfPreviouslyDone (BlastingManager manager, ORF orf, int toProcess) {}
	public void orfRebounded (BlastingManager manager, ORF orf, int start, int stop) {}
	public void orfReboundedAndRename (BlastingManager manager, ORF orf, int start, int stop) {}
	public void orfReincluded(BlastingManager manager, ORF orf) {}
}
