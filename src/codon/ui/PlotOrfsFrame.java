package codon.ui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Desktop;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.List;

import javax.swing.ButtonGroup;
import javax.swing.ImageIcon;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JDialog;
import javax.swing.JEditorPane;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JScrollBar;
import javax.swing.JSeparator;
import javax.swing.JSplitPane;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;
import javax.swing.Timer;
import javax.swing.border.CompoundBorder;
import javax.swing.border.EmptyBorder;
import javax.swing.border.LineBorder;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;
import javax.swing.filechooser.FileNameExtensionFilter;

import codon.Starter;
import codon.blast.BlastingMonitor;
import codon.blast.BlastingOptionDialog;
import codon.data.CodonConfiguration;
import codon.data.CodonConfigurationListener;
import codon.data.GeneFeature;
import codon.data.ORF;
import codon.recognizer.OrfRecognizer;
import codon.ui.FilterBox.UserFilter;
import codon.ui.component.JSelectorMenuItem;

public class PlotOrfsFrame extends JFrame {
	private static final long serialVersionUID = 1L;

	public JMenuItem menuLoadFromCodonProject = new JMenuItem ("Load Codon Project");
	public JMenuItem menuLoadOrfsFromEmbl = new JMenuItem ("Load orfs from .embl / .gb");
	public JMenuItem menuStartFromFasta = new JMenuItem ("Start from .fasta");
	public JMenuItem menuFileSaveCodonProject = new JMenuItem ("Save Codon Project");
	public JMenuItem menuFileSaveAsCodonProject = new JMenuItem ("Save Codon Project as...");
	public JMenuItem menuFileImportOrfsFromProject = new JMenuItem ("Import blasted orf from other Codon project");
	public File currentSaveFile = null;
	public JMenuItem menuFileExport = new JMenuItem ("Export as .embl");
	public JMenuItem menuFileExportCDS = new JMenuItem ("Export ORF as .fasta");
	public JMenuItem menuFileExportCDSAmino = new JMenuItem ("Export Amino ORF as .fasta");
	public JMenuItem menuFileLoadFeatureGene = new JMenuItem ("Load gene feature from .embl");
	public JMenuItem menuFileLoadExit = new JMenuItem ("Exit");
	
	public JCheckBoxMenuItem menuViewHideRemovedOrf = new JCheckBoxMenuItem ("Hide removed orfs");
	
	public JCheckBoxMenuItem menuViewStartStops = new JCheckBoxMenuItem ("Show Start/Stops");
	public JCheckBoxMenuItem menuViewFeatureGene = new JCheckBoxMenuItem ("Show existing gene feature from Embl/GenBank");
	public JCheckBoxMenuItem menuViewValina = new JCheckBoxMenuItem ("Show Valines/Leucines");
	public JCheckBoxMenuItem menuViewQueues = new JCheckBoxMenuItem ("Show queues");
	public JCheckBoxMenuItem menuViewIntergene = new JCheckBoxMenuItem ("Show integene area");
	public JCheckBoxMenuItem menuViewGC = new JCheckBoxMenuItem ("Show GC");
	public JMenuItem menuViewGoTo = new JMenuItem ("Go to ...");
	public JMenuItem menuViewMark = new JMenuItem ("Mark");
	public JRadioButtonMenuItem menuViewBaseSequence = new JRadioButtonMenuItem ("View Base sequence");
	public JRadioButtonMenuItem menuViewAcidoSequence = new JRadioButtonMenuItem ("View Code sequence");
	public JMenuItem menuViewReport = new JMenuItem ("Report");
	
	public JMenuItem menuActionOptimizedBlast = new JMenuItem ("Optimized Blast with uniprot");
	public JMenuItem menuActionFullBlast = new JMenuItem ("Full Blast with uniprot");
	
	public JMenuItem menuActionPrintQueues = new JMenuItem ("TODO - Mark hypotetic broken queues");
//	public JMenuItem menuActionSaveAlteredORF = new JMenuItem ("Save altered orfs");
	
	public JMenuItem menuActionAddRNAOrfs = new JMenuItem ("Annotate RNA");
//	public JMenuItem menuActionAddTRNAOrfs = new JMenuItem ("Match tRNA ORFs");
	
//	public JMenuItem menuFilterRestart = new JMenuItem ("Restart");
	
	public JMenuItem menuFilterReboundStart = new JMenuItem ("Resize using alternative start M/V ");
	public JMenuItem menuFilterSelectBetterFit = new JMenuItem ("Select entries that better fit");
	public JMenuItem menuFilterSelectFulFill = new JMenuItem ("Select entries to fulfill intergenic regions");
//	public JMenuItem menuFilterClear = new JMenuItem ("ORFs: clear filter");
	public JMenuItem menuFilterRemoveTinyUnsignificant = new JMenuItem ("Remove tiny unsignificant orf");
	public JMenuItem menuFilterSuperpositions = new JMenuItem ("Remove overlaps");
	public JSelectorMenuItem menuFilterLowAccuracy = new JSelectorMenuItem ("Remove low accuracy", 80);
	public JMenu menuFilterUserFilter = new JMenu ("User filters");
	public JMenuItem menuFilterCreateFilter = new JMenuItem ("Create filter");
	
	public JMenu menuSelection = new JMenu ("Selection");
	public JCheckBoxMenuItem menuSelectionForceRemoved = new JCheckBoxMenuItem ("Force as removed");
	public JCheckBoxMenuItem menuSelectionForceIncluded = new JCheckBoxMenuItem ("Force as included");
	public JMenuItem menuSelectionForceReboundLeft = new JMenuItem ("Move start to the  left");
	public JMenuItem menuSelectionForceReboundRight = new JMenuItem ("Move start to the  right");
	public JMenuItem menuSelectionAutoRebound = new JMenuItem ("Resize using alternative start M/V");
	public JMenuItem menuSelectionSelectFit = new JMenuItem ("Select Hit entry that better fit");
	public JMenuItem menuSelectionBlast = new JMenuItem ("Blast");
	
	public JMenuItem menuOptionParameters = new JMenuItem ("Parameters");
	public JCheckBoxMenuItem menuOptionTraceJUniprotAPI = new JCheckBoxMenuItem ("Trace UNIPROT API error");
	
	public JMenuItem menuAboutAbout = new JMenuItem ("About");
		
	public PlotOrfsPanel plotOrfsPanel;
	public JPanel statusBar;
	JScrollBar scroll = new JScrollBar (JScrollBar.VERTICAL);
	
	JLabel labThread;
	
	public TabRight tabRight = null;
	
	public OrfRecognizer reco = null;
	
	public boolean changed = true;
	
	public static PlotOrfsFrame instance = null;
	
	ActionListener repaintActionListener = new ActionListener () {
		public void actionPerformed (ActionEvent e) {
			plotOrfsPanel.recompute ();
		}
	};
	
	JFileChooser chooser = new JFileChooser ();
		
	public PlotOrfsFrame (OrfRecognizer _reco, boolean startBlastingQueue, boolean loadEveryOrfsDetailsAtStart, boolean loadRNAAtStart) {
		this.reco = _reco;
		instance = this;
		setTitle ("CODON - " + reco.dataDestDirS + " * orfs:" + reco.candidateOrfCount ());
		
		setIconImage (Starter.icon);
		
		JMenuBar bar = new JMenuBar ();
		JMenu menuFile = new JMenu ("File");
		JMenu menuView = new JMenu ("View");
		JMenu menuAction = new JMenu ("Action");
		JMenu menuFilter = new JMenu ("Curation");
		JMenu menuOption = new JMenu ("Options");
		JMenu menuAbout = new JMenu ("About");
		bar.add (menuFile);
		bar.add (menuView);
		bar.add (menuAction);
		bar.add (menuFilter);
		bar.add (menuSelection);
		bar.add (menuOption);
		bar.add (menuAbout);
		setJMenuBar (bar);
		
		menuFile.add (menuLoadFromCodonProject);
		menuLoadFromCodonProject.setAccelerator (KeyStroke.getKeyStroke ("control O"));
		menuFile.add (menuLoadOrfsFromEmbl);
		menuFile.add (menuStartFromFasta);
		menuFile.add (new JSeparator ());
		menuFile.add (menuFileSaveCodonProject);
		menuFileSaveCodonProject.setAccelerator (KeyStroke.getKeyStroke ("control S"));
		menuFile.add (menuFileSaveAsCodonProject);
		menuFileSaveAsCodonProject.setAccelerator (KeyStroke.getKeyStroke ("control shift S"));
		menuFile.add (new JSeparator ());
		menuFile.add (menuFileImportOrfsFromProject);
		menuFile.add (new JSeparator ());
		menuFile.add (menuFileExport);
		menuFile.add (menuFileExportCDS);
		menuFile.add (menuFileExportCDSAmino);
//		menuFile.add (menuFileLoadFeatureGene);
		menuFile.add (menuFileLoadExit);
//		menuFileExport.setAccelerator (KeyStroke.getKeyStroke ("control E"));
		
		menuView.add (menuViewStartStops);
		menuViewStartStops.setAccelerator (KeyStroke.getKeyStroke ("S"));
		menuView.add (menuViewValina);
		menuViewValina.setAccelerator (KeyStroke.getKeyStroke ("V"));
		menuView.add (menuViewHideRemovedOrf);
		menuViewHideRemovedOrf.setAccelerator (KeyStroke.getKeyStroke ("control R"));
		menuViewHideRemovedOrf.setSelected (true);
		menuView.add (menuViewFeatureGene);
		menuViewFeatureGene.setSelected (false);
		menuView.add (new JSeparator ());
//		menuView.add (menuViewQueues);
//		menuView.add (menuViewIntergene);
		menuViewFeatureGene.setSelected (false);
		menuViewIntergene.setSelected (true);
		menuViewQueues.setAccelerator (KeyStroke.getKeyStroke ("control Q"));
		menuViewQueues.setSelected (true);
//		menuView.add (menuViewGC);
		menuView.add (menuViewBaseSequence);
		menuView.add (menuViewAcidoSequence);
		menuView.add (new JSeparator ());
		menuView.add (menuViewGoTo);
		menuViewGoTo.setAccelerator (KeyStroke.getKeyStroke ("control G"));
		menuView.add (menuViewMark);
		menuViewMark.setAccelerator (KeyStroke.getKeyStroke ("control M"));
		menuView.add (new JSeparator ());
		menuView.add (menuViewReport);
		
		menuAction.add (menuActionOptimizedBlast);
		menuAction.add (menuActionFullBlast);
//		menuAction.add (menuActionUpgradeWithSwiss);
		menuAction.add (new JSeparator ());
		menuAction.add (menuActionAddRNAOrfs);
//		menuAction.add (menuActionAddTRNAOrfs);
//		menuAction.add (new JSeparator ());
//		menuAction.add (menuActionPrintQueues);
		
//		menuFilter.add (menuFilterRestart);
		menuFilter.add (menuFilterReboundStart);
		menuFilterReboundStart.setAccelerator (KeyStroke.getKeyStroke ("control alt A"));
		menuFilter.add (menuFilterSelectBetterFit);
		menuFilterSelectBetterFit.setAccelerator (KeyStroke.getKeyStroke ("control alt F"));
		menuFilter.add (menuFilterSelectFulFill);
		menuFilter.add (new JSeparator ());
//		menuFilter.add (menuFilterClear);
//		menuFilter.add (menuFilterRemoveTinyUnsignificant);
		menuFilter.add (menuFilterLowAccuracy);
		menuFilter.add (menuFilterSuperpositions);
		menuFilter.add (new JSeparator ());
		menuFilter.add (menuFilterUserFilter);
		buildFilterList ();
		
		menuSelection.add (menuSelectionForceRemoved);
		menuSelectionForceRemoved.setAccelerator (KeyStroke.getKeyStroke ("DELETE"));
		menuSelection.add (menuSelectionForceIncluded);
		menuSelectionForceIncluded.setAccelerator (KeyStroke.getKeyStroke ("INSERT"));
		menuSelection.add (new JSeparator ());
		menuSelection.add (menuSelectionForceReboundLeft);
		menuSelectionForceReboundLeft.setAccelerator (KeyStroke.getKeyStroke ("control LEFT"));
		menuSelection.add (menuSelectionForceReboundRight);
		menuSelectionForceReboundRight.setAccelerator (KeyStroke.getKeyStroke ("control RIGHT"));
		menuSelection.add (menuSelectionAutoRebound);
		menuSelection.add (new JSeparator ());
		menuSelectionAutoRebound.setAccelerator (KeyStroke.getKeyStroke ("control A"));
		menuSelection.add (menuSelectionSelectFit);
		menuSelectionSelectFit.setAccelerator (KeyStroke.getKeyStroke ("control F"));
		menuSelection.add (new JSeparator ());
		menuSelection.add (menuSelectionBlast);
		menuSelectionBlast.setAccelerator (KeyStroke.getKeyStroke ("control B"));
		
		menuOption.add (menuOptionParameters);
		menuOption.add (new JSeparator ());
		menuOption.add (menuOptionTraceJUniprotAPI);
		
		ButtonGroup radioSequenceGroup = new ButtonGroup ();
		radioSequenceGroup.add (menuViewBaseSequence);
		menuViewAcidoSequence.setSelected (true);
		radioSequenceGroup.add (menuViewAcidoSequence);
		
		plotOrfsPanel = new PlotOrfsPanel (reco, this);
		
		menuViewStartStops.addActionListener (repaintActionListener);
		menuViewValina.addActionListener (repaintActionListener);
		menuViewBaseSequence.addActionListener (repaintActionListener);
		menuViewAcidoSequence.addActionListener (repaintActionListener);
		menuViewHideRemovedOrf.addActionListener (repaintActionListener);
		menuViewQueues.addActionListener (repaintActionListener);
		menuViewGC.addActionListener (repaintActionListener);
		menuViewIntergene.addActionListener (repaintActionListener);
		menuViewFeatureGene.addActionListener (repaintActionListener);
		
		menuAbout.add (menuAboutAbout);

		menuFileExport.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				exportEmbl ();
			}
		});
		menuFileExportCDS.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent e) {
				try {
					chooser.setSelectedFile (currentSaveFile);
					FileNameExtensionFilter filter = new FileNameExtensionFilter ("Fasta", "fasta", "fa");
					chooser.setFileFilter (filter);
					if (chooser.showSaveDialog (PlotOrfsFrame.this) == JFileChooser.APPROVE_OPTION) {
						BufferedWriter writer = new BufferedWriter (new FileWriter (chooser.getSelectedFile ()));
						for (int i = 0; i < reco.orfsPorFrame.length; i++) {
							for (ORF orf1: reco.orfsPorFrame [i]) {
								if (orf1.isRemoved () || orf1.product != null && orf1.product.entryId != null && (orf1.product.entryId.equals ("RNA") || orf1.product.entryId.equals ("tRNA"))) continue; 
								writer.write (">" + orf1.id);
								writer.newLine ();
								writer.write (orf1.getSequence ());
								writer.newLine ();
							}
						}
						writer.close ();
					}
				}
				catch (IOException ex) {
					ex.printStackTrace ();
				}
			}
		});
		menuFileExportCDSAmino.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent e) {
				try {
					chooser.setSelectedFile (currentSaveFile);
					FileNameExtensionFilter filter = new FileNameExtensionFilter ("Fasta", "fasta", "fa");
					chooser.setFileFilter (filter);
					if (chooser.showSaveDialog (PlotOrfsFrame.this) == JFileChooser.APPROVE_OPTION) {
						BufferedWriter writer = new BufferedWriter (new FileWriter (chooser.getSelectedFile ()));
						for (int i = 0; i < reco.orfsPorFrame.length; i++) {
							for (ORF orf1: reco.orfsPorFrame [i]) {
								if (orf1.isRemoved () || orf1.product != null && orf1.product.entryId != null && (orf1.product.entryId.equals ("RNA") || orf1.product.entryId.equals ("tRNA"))) continue; 
								writer.write (">" + orf1.id);
								writer.newLine ();
								writer.write (orf1.getAminoSequence ());
								writer.newLine ();
							}
						}
						writer.close ();
					}
				}
				catch (IOException ex) {
					ex.printStackTrace ();
				}
			}
		});
		menuFileSaveCodonProject.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent e) {
				saveProject (false);
			}
		});
		menuFileImportOrfsFromProject.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				importOrf ();
			}
		});
		menuFileSaveAsCodonProject.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent e) {
				saveAsProject ();
			}
		});
		menuLoadFromCodonProject.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent e) {
				loadProject ();
			}
		});
		menuLoadOrfsFromEmbl.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent e) {
				loadFromEmbl ();
			}
		});
		menuStartFromFasta.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				startFromFasta ();
			}
		});
		menuFileLoadFeatureGene.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				chooser.setCurrentDirectory (new File (CodonConfiguration.workspace));
				FileNameExtensionFilter filter = new FileNameExtensionFilter ("Embl", "embl");
				chooser.setFileFilter (filter);
				if (chooser.showOpenDialog (PlotOrfsFrame.this) == JFileChooser.APPROVE_OPTION) {
					loadGeneFeature (chooser.getSelectedFile ());
				}
			}
		});
		menuFileLoadExit.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				quit ();
			}
		});
		
		menuViewGoTo.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				String v = null;
				if ((v = JOptionPane.showInputDialog (PlotOrfsFrame.this, "Enter the position sequence to display")) != null) {
					try {
						int pos = Integer.parseInt (v);
						plotOrfsPanel.setSelectedPosition (pos, true);
						return;
					} catch (Exception e) {
					}
				}
			}
		});
		menuViewReport.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				new ReportDialog (PlotOrfsFrame.this);
			}
		});
		
		menuActionOptimizedBlast.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent e) {
				new BlastingOptionDialog (PlotOrfsFrame.this, reco, false);
			}
		});
		menuActionFullBlast.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent e) {
				new BlastingOptionDialog (PlotOrfsFrame.this, reco, true);
			}
		});

		menuActionPrintQueues.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent e) {
				try {
					BufferedWriter writer = new BufferedWriter (new FileWriter (CodonConfiguration.workspace + reco.dataDestDirS + "/queues.txt"));
					for (List <ORF> orfs: reco.orfsPorFrame) {
						for (ORF orf: orfs) {
							if (orf.isRemoved ()) continue;
							orf.printTail (writer);
						}
					}
					writer.close ();
				} catch (IOException e1) {
					e1.printStackTrace ();
				}
			}
		});
		menuFilterReboundStart.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				reco.reboundStarts ();
				plotOrfsPanel.recompute ();
				changed = true;
				updateFilter ();
			}
		});
		menuFilterSelectBetterFit.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				reco.selectBetterEntriesThatFit ();
				plotOrfsPanel.recompute ();
				changed = true;
				updateFilter ();
			}
		});
		menuFilterSelectFulFill.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				reco.selectBetterEntriesFulfillIntergenic (tabRight.overlapPane.getTolerance ());
				plotOrfsPanel.recompute ();
				changed = true;
				updateFilter ();
			}
		});
		menuViewMark.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				addMark ();
			}
		});
//		menuActionSaveAlteredORF.addActionListener (new ActionListener () {
//			public void actionPerformed (ActionEvent arg0) {
//				int renamed = 0;
//				for (int i = 0; i < 6; i++) {
//					for (ORF orf: reco.orfsPorFrame [i]) {
//						if (orf.isRemoved ()) continue;
//						String lg = "" + orf.lenght;
//						for (int k = lg.length (); k < 6; k ++) lg = '0' + lg; 
//						String novoName = "F" + i + "-" + lg + "-" + orf.start + "-" + orf.stop;
//						if (!novoName.equals (orf.id)) {
//							renamed ++;	
//							orf.saveAlterations (novoName);
//						}
//					}
//				}
//				System.err.println (renamed + " files renamed");
//			}
//		});
//		menuFilterClear.addActionListener (new ActionListener () {
//			public void actionPerformed (ActionEvent e) {
//				reco.clearRemoved ();
//				menuFilterRemoveTinyUnsignificant.setEnabled (true);
//				menuFilterSuperpositions.setEnabled (true);
//				menuFilterLowAccuracy.setEnabled (true);
//				updateFilter ();
//			}
//		});

//		menuFilterRemoveTinyUnsignificant.addActionListener (new ActionListener () {
//			public void actionPerformed (ActionEvent e) {
//				reco.removeTinyUnspecifyingSequence ();
//				changed = true;
//				updateFilter ();
//			}
//		});
		menuFilterSuperpositions.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent e) {
				reco.removeOverlaps (tabRight.overlapPane.getTolerance ());
				changed = true;
				updateFilter ();
			}
		});
		menuFilterLowAccuracy.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent e) {
				reco.removeLowAccuracy (menuFilterLowAccuracy.getValue ());
//				menuFilterLowAccuracy.setEnabled (false);
				changed = true;
				updateFilter ();
			}
		});
		menuFilterCreateFilter.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent e) {
				new FilterBox (PlotOrfsFrame.this, reco);
			}
		});
//		menuFilterRestart.addActionListener (new ActionListener () {
//			public void actionPerformed (ActionEvent e) {
//				menuFilterRemoveTinyUnsignificant.setEnabled (true);
//				menuFilterSuperpositions.setEnabled (true);
//				menuFilterLowAccuracy.setEnabled (true);
//				if (!reco.loadedFromProject) {
//					reco.clear ();
//					reco.load (reco.dataSrcFileS);
//					loadORFs ();
//				}
//				else {
//					reco.clearRemoved ();
//				}
//				updateFilter ();
//			}
//		});
		menuActionAddRNAOrfs.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent e) {
				reco.annotateRNAOrf (true);
				changed = true;
				updateFilter ();
			}
		});
//		menuActionAddTRNAOrfs.addActionListener (new ActionListener () {
//			public void actionPerformed (ActionEvent e) {
//				reco.addTRNAOrf (true);
//				changed = true;
//				updateFilter ();
//			}
//		});
		menuSelectionForceRemoved.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				boolean excluing = menuSelectionForceRemoved.isSelected ();
				plotOrfsPanel.getSelected ().orf.setExcluded (excluing);
				if (excluing) {
					plotOrfsPanel.setSelected (null);
					menuSelectionForceIncluded.setSelected (false);
				}
				changed = true;
				updateFilter ();
			}
		});
		menuSelectionForceIncluded.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				boolean included = menuSelectionForceIncluded.isSelected ();
				plotOrfsPanel.getSelected ().orf.setForced (included);
				if (included) {
					menuSelectionForceRemoved.setSelected (false);
				}
				changed = true;
				updateFilter ();
			}
		});
		menuSelectionForceReboundLeft.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				plotOrfsPanel.getSelected ().orf.reboundLeft (reco.startStopPorFrame [plotOrfsPanel.getSelected ().orf.frame]);
				changed = true;
				plotOrfsPanel.repaint ();
//				updateFilter ();
			}
		});
		menuSelectionForceReboundRight.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				plotOrfsPanel.getSelected ().orf.reboundRight (reco.startStopPorFrame [plotOrfsPanel.getSelected ().orf.frame]);
				changed = true;
				plotOrfsPanel.repaint ();
//				updateFilter ();
			}
		});
		menuSelectionAutoRebound.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				plotOrfsPanel.getSelected ().orf.autoRebound ();
				changed = true;
				updateFilter ();
			}
		});
		menuSelectionSelectFit.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				plotOrfsPanel.getSelected ().orf.selectBetterEntryThatFit ();
				changed = true;
				updateFilter ();
			}
		});
		menuSelectionBlast.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				if (plotOrfsPanel.getSelected ().orf != null) tabRight.addBlast (plotOrfsPanel.getSelected ().orf);
			}
		});
		menuOptionParameters.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				new CodonConfiguration (PlotOrfsFrame.this);
			}
		});
		menuOptionTraceJUniprotAPI.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				BlastingMonitor.blastErrors = menuOptionTraceJUniprotAPI.isSelected ();
			}
		});
		menuAboutAbout.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent e) {
//				String mess = "<HTML><BODY> "
//						+ "This project is sustained by the following volunteer associated members (by name): "
//						+ "<UL> "
//						+ "<LI>Adonney Allan de Oliveira Veras (Phd in Bioinformatics, UFPA - <BR>NDAE/PPCA) </LI>"
//						+ "<LI>Bruno Merlin (Phd in Computation/HCI, UFPA - NDAE/PPCA)</LI>"
//						+ "<LI>Jorianne Thyeska Castro Alves (Phd in Genetic and Molecular Biology, UEPA)</LI>"
//						+ "<Li>M�nica Silva de Oliveira (Master in Computation, UFPA - NDAE/PPCA)</LI>"
//						+ "<LI>Pablo Henrique Caracciolo Gomes de S� (Phd in Bioinformatics, UFRA)</LI>"
//						+ "</UL>"
//						+ "<BR>No member is earning any money keeping this project alive (and never will)."
//						+ " However, we <BR>need your help to offer you a support, to correct the bugs, "
//						+ "to implement new fancy features <BR>(and to keep alive our poor public infrastructure). "
//						+ "Then, if this project is helping you in<BR> anything, please, help us! "
//						+ "<BR><BR><center><img src=\"file:data/img/cat.png\"/></center> <br>"
//						+ "<BR> PS: If you want to congratulate us, you can send e email direct to me."
//						+ "However, <BR> for any question, bug report or new requirement, "
//						+ "send an email to <A href=\"codon.project@gmail\"> codon.project@gmail.com </A><BR> (I�m currently enjoying holidays!)." + 
//						"</BODY></HTML>";
				String mess = "<HTML><BODY> "
						+ "This project is sustained by the following volunteer associated members (by name): "
						+ "<UL> "
						+ "<LI>Adonney Allan de Oliveira Veras (Phd in Bioinformatics, UFPA - <BR>NDAE/PPCA) </LI>"
						+ "<LI>Bruno Merlin (Phd in Computation/HCI, UFPA - NDAE/PPCA)</LI>"
						+ "<LI>Jorianne Thyeska Castro Alves (Phd in Genetic and Molecular Biology, UEPA)</LI>"
						+ "<Li>M�nica Silva de Oliveira (Master in Computation, UFPA - NDAE/PPCA)</LI>"
						+ "<LI>Pablo Henrique Caracciolo Gomes de S� (Phd in Bioinformatics, UFRA)</LI>"
						+ "</UL>"
						+ "<BR> For any question, bug report or new requirement, "
						+ "send an email to <A href=\"codon.project@gmail\"> codon.project@gmail.com </A><BR> " + 
						"</BODY></HTML>";
//				JOptionPane.showMessageDialog (PlotOrfsFrame.this, mess, "Project Team", JOptionPane.INFORMATION_MESSAGE, new ImageIcon ("data/img/reloadA.png"));
				JEditorPane ep;
				ep = new JEditorPane ("text/html", mess);
				ep.setEditable (false);
				ep.addHyperlinkListener (new HyperlinkListener () {
					public void hyperlinkUpdate (HyperlinkEvent e) {
						if (e.getEventType() == HyperlinkEvent.EventType.ACTIVATED) {
							try {
								Desktop.getDesktop ().browse (e.getURL ().toURI ());
							} catch (IOException | URISyntaxException ex) {
							}
						}
			        }
			    });
//				ep.setPreferredSize (new Dimension (ep.getPreferredSize ().width, ep.getPreferredSize ().height + 115));
			    Object [] buttons = {"Ok"};
			    JOptionPane.showOptionDialog (PlotOrfsFrame.this, ep, "CODON Project Team - version 1.2 Beta 4", JOptionPane.DEFAULT_OPTION, JOptionPane.QUESTION_MESSAGE, new ImageIcon ("data/img/Logo1.png"), buttons, buttons[0]);
			}
		});
		
		statusBar = new JPanel (new FlowLayout (FlowLayout.LEFT));
	    statusBar.setBorder (new CompoundBorder (new LineBorder (Color.DARK_GRAY), new EmptyBorder (4, 4, 4, 4)));
	    
	    JSplitPane splitPane = new JSplitPane (JSplitPane.HORIZONTAL_SPLIT);
	    
	    JPanel plotMain = new JPanel (new BorderLayout ());
	    plotMain.add (plotOrfsPanel);
	    plotMain.add (scroll, BorderLayout.EAST);
	    plotMain.add (new SequenceEditor (this), BorderLayout.SOUTH);
	    splitPane.add (plotMain);
	    tabRight = new TabRight (this);
	    splitPane.add (tabRight);
	    
	    getContentPane ().add (splitPane);
	    getContentPane ().add (statusBar, BorderLayout.SOUTH);
	    labThread = new JLabel ();
	    updateThread ();
	    statusBar.add (labThread);

	    scroll.setMinimum (0);
	    scroll.setValue (0);
	    scroll.setMaximum (reco.length ());
	    scroll.addAdjustmentListener (new AdjustmentListener () {
			public void adjustmentValueChanged (AdjustmentEvent ev) {
				plotOrfsPanel.setSelectedPosition (ev.getValue (), false);
			}
		});
	    configureScroll ();
	    plotOrfsPanel.addComponentListener (new ComponentAdapter () {
	    	public void componentResized (ComponentEvent e) {
	    		configureScroll ();
	    		splitPane.setDividerLocation (getWidth () - 250);
	    	}
		});
//		setDefaultCloseOperation (JFrame.EXIT_ON_CLOSE);
	    setDefaultCloseOperation (JFrame.DO_NOTHING_ON_CLOSE);
	    addWindowListener (new WindowAdapter () {
	    	public void windowClosing (WindowEvent e) {
	    		quit ();
	    	}
		});
		setSize (1200, 900);
//		if (startBlastingQueue) blastingManager.loadBlastingProcesses ();
		setExtendedState (JFrame.MAXIMIZED_BOTH);
		if (loadEveryOrfsDetailsAtStart) loadORFs (loadRNAAtStart);
		chooser.setCurrentDirectory (new File (CodonConfiguration.workspace));
		
		Timer monitor = new Timer (10000, new ActionListener () {
			public void actionPerformed (ActionEvent e) {
				updateThread ();
			}
		});
		monitor.setRepeats (true);
		monitor.start ();
		
		CodonConfiguration.addCodonConfigurationListener (new CodonConfigurationListener () {
			public void minOrfLengthChanged (int minOrfLength) {
				for (ArrayList <ORF> orfs: reco.orfsPorFrame) for (ORF orf: orfs) {
					if (orf.length < CodonConfiguration.min_orf_length) {
						if (!orf.isBlasted () 
								|| orf.getAccuracy () < CodonConfiguration.accuracy_level_to_remove_blasted_orfs_during_blast
								|| orf.isUncharacterized ()) {
							orf.setRemoved (true);
						}
					}
				}
			}
		});
		
		setVisible (true);
	}
	
	protected void quit () {
		if (changed) {
			int v = JOptionPane.showConfirmDialog (this, "Do you want to save the alteration as a CODON project?");
			if (JOptionPane.OK_OPTION == v) {
				if (saveAsProject ()) System.exit (0);
			}
			else if (JOptionPane.NO_OPTION == v) {
				System.exit (0);
			}
		}
		else {
			System.exit (0);
		}
		
	}

	protected void importOrf () {
		OrfRecognizer reference = new OrfRecognizer ();
		chooser.setCurrentDirectory (new File (CodonConfiguration.workspace));
		FileNameExtensionFilter filter = new FileNameExtensionFilter ("Codon", "codon");
		chooser.setFileFilter (filter);
		if (chooser.showOpenDialog (PlotOrfsFrame.this) == JFileChooser.APPROVE_OPTION) {
			reference.loadFromCodonProject (chooser.getSelectedFile ());
			new ImportORFSplash (this, reco, reference);
		}
		changed = true;
		updateFilter ();
	}

	public void updateFilter () {
		reco.computeOrfFreeAndOrfSuperpositionSpans (tabRight.overlapPane.getTolerance ()); 
		plotOrfsPanel.recompute ();
		tabRight.organismPane.updateFamillies ();
		updateTitle ();
	}
	
	public void updateTitle () {
		reco.computeOrfFreeAndOrfSuperpositionSpans (tabRight.overlapPane.getTolerance ()); 
		plotOrfsPanel.recompute ();
		try {
			setTitle ("CODON - ORF src: " + reco.dataDestDirS + "      " + (currentSaveFile != null ? currentSaveFile.getCanonicalPath () : "") + (changed? "*" : "") + " orfs: " + reco.candidateOrfCount ());
		} catch (IOException e) {
			e.printStackTrace ();
		}
	}
	
	protected boolean saveAsProject () {
		return saveProject (true);
	}
	
	protected boolean saveProject (boolean forceNew) {
		if (currentSaveFile == null || forceNew) {
			currentSaveFile = new File (CodonConfiguration.workspace + reco.dataDestDirS + "/auto.codon");
			chooser.setCurrentDirectory (new File (CodonConfiguration.workspace));
			chooser.setSelectedFile (currentSaveFile);
			FileNameExtensionFilter filter = new FileNameExtensionFilter ("Codon", "codon");
			chooser.setFileFilter (filter);
			if (chooser.showSaveDialog (PlotOrfsFrame.this) == JFileChooser.APPROVE_OPTION) {
				currentSaveFile = chooser.getSelectedFile ();
			}
			else return false;
		}
		reco.saveAsCodonProject (currentSaveFile);
		changed = false;
		updateTitle ();
		return true;
	}
	
	protected void loadProject () {
		chooser.setCurrentDirectory (new File (CodonConfiguration.workspace));
		FileNameExtensionFilter filter = new FileNameExtensionFilter ("Codon", "codon");
		chooser.setFileFilter (filter);
		if (chooser.showOpenDialog (PlotOrfsFrame.this) == JFileChooser.APPROVE_OPTION) {
			currentSaveFile = chooser.getSelectedFile ();
			reco = new OrfRecognizer ();
			reco.loadFromCodonProject (currentSaveFile);
			plotOrfsPanel.reco = reco;
			changed = false;
			updateFilter ();
		}
	}

	protected void startFromFasta () {
		FileNameExtensionFilter filter = new FileNameExtensionFilter ("Fasta", "fasta", "fa");
		chooser.setFileFilter (filter);
		if (chooser.showOpenDialog (PlotOrfsFrame.this) == JFileChooser.APPROVE_OPTION) {
			File f= chooser.getSelectedFile ();
			currentSaveFile = null;
			reco = new OrfRecognizer ();
			reco.load (f.getAbsolutePath ());
			plotOrfsPanel.reco = reco;
			changed = true;
			updateFilter ();
		}
		loadORFs (true);
	}
	
	protected void loadFromEmbl () {
		FileNameExtensionFilter filter = new FileNameExtensionFilter ("Embl & GenBank", "embl", "gb");
		chooser.setFileFilter (filter);
		if (chooser.showOpenDialog (PlotOrfsFrame.this) == JFileChooser.APPROVE_OPTION) {
			reco = new OrfRecognizer ();
			File file = chooser.getSelectedFile ();
			currentSaveFile = null;
			reco.loadFromEmbl (file);
			plotOrfsPanel.reco = reco;
			changed = true;
			updateFilter ();
			sendRNAS ();
		}
		loadORFs (false);
	}
	
//	class ErrorControlerReader extends Thread {
//		public void run () {
//			InputStream inP = proc.getErrorStream ();
//			BufferedReader reader = new BufferedReader (new InputStreamReader (inP));
//			String line;
//			try {
//				while ((line = reader.readLine ()) != null) {
//					System.err.println ("ERROR: " + line);
//				}
//			}
//			catch (IOException io) {
//			}
//		}
//	}
//	Process proc;
	public void sendRNAS () {
		reco.upgradeRNAS ();
	}
	
	protected void exportEmbl () {
		chooser.setCurrentDirectory (new File (CodonConfiguration.workspace + reco.dataDestDirS));
		FileNameExtensionFilter filter = new FileNameExtensionFilter ("Embl", "embl");
		chooser.setFileFilter (filter);
		chooser.setSelectedFile (new File (CodonConfiguration.workspace + reco.dataDestDirS + "/" + reco.dataDestDirS + ".embl"));
		if (chooser.showSaveDialog (PlotOrfsFrame.this) != JFileChooser.APPROVE_OPTION) {
			return;
		}
		reco.exportAsEmbl (chooser.getSelectedFile ());
	}

	public void configureScroll () {
		scroll.setUnitIncrement (plotOrfsPanel.getWidth ());
		scroll.setBlockIncrement ((int) (plotOrfsPanel.getWidth () *  plotOrfsPanel.getHeight () / (7. * plotOrfsPanel.scaleY)));
	}

	public void addMark () {
		String comment = JOptionPane.showInputDialog (this, "Add a comment:");
		if (plotOrfsPanel.getSelected () != null) {
			ORF orf = plotOrfsPanel.getSelected ().orf;
			comment = orf.start + ": " + orf.product.productId + " / " + orf.product.entryId + (comment == null ? "" : " - " + comment);
		}
		else {
			comment = plotOrfsPanel.getVisibleIndex () + ": " + (comment == null ? "" : " - " + comment);
		}
		tabRight.addMark (comment);
	}
	
	class LoadOrfDialog extends JDialog {
		private static final long serialVersionUID = 1L;
		int total = 0;
		int done = 0;
		int orfWithoutEntries = 0;
		int orfUnblasted = 0;
		int orfWithoutMatchingSequence = 0;
		JLabel lLoad = new JLabel ("Loading orf 0/" + total);
		JProgressBar progress;
		public LoadOrfDialog (boolean loadRna) {
			for (List <ORF> orfs: reco.orfsPorFrame) {
				for (ORF orf: orfs) {
					if (!orf.isRemoved ()) total ++;
				}
			}
			if (total == 0) return;
			lLoad.setText ("Loading orf " + 1 + " / " + total);
			getContentPane ().add (lLoad);
			
			progress = new JProgressBar (0, total ++);
			getContentPane ().add (progress, BorderLayout.SOUTH);
			getContentPane ().add (new JLabel (new ImageIcon (Starter.iconBordrered)), BorderLayout.NORTH);
			getContentPane ().add (new JPanel (), BorderLayout.WEST);
			getContentPane ().add (new JPanel (), BorderLayout.EAST);
			setUndecorated (true);
			setLocationRelativeTo (PlotOrfsFrame.this);
			setSize (200, 140);
			setModal (true);
			
			Thread loader = new Thread () {
				public void run () {
					for (List <ORF> orfs: reco.orfsPorFrame) {
						for (ORF orf: orfs) {
							if (!orf.isRemoved ()) {
//								System.err.println ("Load " + orf.id);
								orf.loadProduct ();
								if (orf.product == null) {
									System.err.println ("Load pb " + orf.id);
									continue;
								}
								if (!orf.product.blasted) orfUnblasted ++;
								if (!orf.product.hasEntries) orfWithoutEntries ++;
								if (!orf.product.hasMatchingSequence ()) orfWithoutMatchingSequence ++;
								if (orf.product.blasted && orf.product.hasEntries) reco.duplicatedEntries.addString (orf.product.entryId, orf);
								if (orf.product.blasted && orf.product.hasEntries && !orf.product.productId.equals ("?")) reco.duplicatedGenes.addString (orf.product.productId, orf);
//								System.err.println(orf.id);
//								if (!orf.removed && orf.product.hasEntries) orf.product.store (reco.dataDestFileS, orf.id);
								done ++;
								SwingUtilities.invokeLater (new Runnable () {
									public void run () {
										lLoad.setText ("Loading orf " + (done + 1) + " / " + total);
										progress.setValue (done);
									}
								});
							}
						}
					}
					if (loadRna) {
						File f = new File (CodonConfiguration.workspace + reco.dataDestDirS + "/rna.txt");
						if (f.exists ()) {
							menuActionAddRNAOrfs.setText ("Re-annotate RNA");
							reco.annotateRNAOrf (false);
						}
//						f = new File (CodonConfiguration.workspace + reco.dataDestDirS + "/trna.txt");
//						if (f.exists ()) {
//							menuActionAddTRNAOrfs.setText ("Re-match tRNA CDS");
//							reco.addTRNAOrf (false);
//						}
					}
					reco.loadPreannotated ();
					setVisible (false);
//					if (orfUnblasted + orfWithoutEntries + orfWithoutMatchingSequence != 0) {
//						JOptionPane.showMessageDialog (PlotOrfsFrame.this, "Unblasted: " + orfUnblasted 
//																		+ "\nOrfs without entries: " + orfWithoutEntries
//																		+ "\nOrf without matching sequence : " + orfWithoutMatchingSequence);
//					}
//					tabRight.organismPane.updateFamillies ();
//					reco.computeOrfFreeAndOrfSuperpositionSpans (tabRight.overlapPane.getTolerance ()); 
//					plotOrfsPanel.recompute ();
					SwingUtilities.invokeLater (new Runnable () {
						public void run  () {
							updateFilter ();
						}
					});
					
				};
			};
			
//			SwingWorker <Object, Object> worker = new SwingWorker <Object, Object> () {
//				protected Object doInBackground () throws Exception {
//					return null;
//				}
//			};
//			worker.execute ();
			loader.start ();
			setDefaultCloseOperation (DISPOSE_ON_CLOSE);
			setVisible (true);
		}
	}
	
	public void loadORFs (boolean loadRna) {
		new LoadOrfDialog (loadRna);
	}
	
	public void loadGeneFeature (File file) {
		try {
			reco.geneFeatures.clear ();
			BufferedReader reader = new BufferedReader (new FileReader (file));
			String line = "";
			while ((line = reader.readLine ()) != null) {
				boolean reverse = false;
				boolean pseudo = false;
				if (line.contains ("FT   gene")) {
					if (line.contains (">")) {
						line = line.replace (">", "");
						pseudo = true;
					}
					if (line.contains ("<")) {
						line = line.replace ("<", "");
						pseudo = true;
					}
					if (line.contains ("complement")) {
						line = line.replace ("complement(", "");
						line = line.replace (")", "");
						reverse = true;
					}
					String start = line.substring (21, line.indexOf ('.'));
					String stop = line.substring (line.indexOf ('.') + 2);
					String gene = reader.readLine ();
					gene = gene.substring (gene.indexOf ('"') + 1, gene.length () - 1);
					reco.geneFeatures.add (new GeneFeature (gene, Integer.parseInt (start) - 1, Integer.parseInt (stop) - 3, reverse, pseudo));
				}
			}
			reader.close ();
		} catch (IOException e) {
			e.printStackTrace ();
		}
		plotOrfsPanel.recompute ();
	}

	public void updateThread () {
		MemoryMXBean memoryMXBean = ManagementFactory.getMemoryMXBean ();
		labThread.setText (Thread.activeCount () + " running threads - Memory usage" + String.format (": %.2f GB", (double) memoryMXBean.getHeapMemoryUsage ().getUsed () / 1073741824));
	}
	
	public void buildFilterList () {
		menuFilterUserFilter.removeAll ();
		menuFilterUserFilter.add (menuFilterCreateFilter);
		menuFilterUserFilter.add (new JSeparator ());
		for (UserFilter f: FilterBox.savedFilters) {
			JMenuItem m = new JMenuItem (f.name);
			m.addActionListener (new ActionListener () {
				public void actionPerformed (ActionEvent e) {
					FilterBox.applyFilter (PlotOrfsFrame.this, f.name);
					updateFilter ();
				}
			});
			menuFilterUserFilter.add (m);
		}
	}
		
}
