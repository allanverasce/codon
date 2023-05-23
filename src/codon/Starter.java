package codon;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

import javax.imageio.ImageIO;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField; 
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.filechooser.FileNameExtensionFilter;

import codon.data.CodonConfiguration;
import codon.recognizer.OrfRecognizer;
import codon.ui.PlotOrfsFrame;

public class Starter extends JDialog {
	private static final long serialVersionUID = 1L;
	
	public static Image icon = null;
	public static Image iconBordrered = null;
	
	static {
		File dataRep = new File ("data/config");
		dataRep.mkdirs ();
		
		try {
			icon = ImageIO.read (new File ("./data/img/Logo1.png"));
			iconBordrered = ImageIO.read (new File ("./data/img/BLogo.png"));
		} catch (IOException e2) {
			e2.printStackTrace ();
		}
	}

	JTextField fFastaFile = new JTextField ("");
	JTextField fCodonFile = new JTextField ("");
	JTextField fEmblFile = new JTextField ("");	
	
	JTextField fWorkspace = new JTextField ("workspace/");
	JButton bOpenWorkspace = new JButton ("...");

	JButton bOpenChooserFasta = new JButton ("...");
	JButton bOpenChooserCodon = new JButton ("...");
	JButton bOpenChooserEmbl = new JButton ("...");
	JButton bStart = new JButton ("Start");
	JFileChooser chooser = new JFileChooser ();
	
	JPanel panelFileChooseFasta = new JPanel ();
	JPanel panelFileChooseCodon = new JPanel ();
	JPanel panelFileChooseEmbl = new JPanel ();
	
	JCheckBox ckStartBlastingQueue = new JCheckBox ("Start blasting queue");
	JCheckBox ckLoadEveryORFDetails = new JCheckBox ("Start loading every orf details");
	
	JRadioButton radioCodon = new JRadioButton ("Codon project");
	JRadioButton radioFasta = new JRadioButton ("Fasta data");
	JRadioButton radioEmbl = new JRadioButton ("Embl/GB file");
	
	public Starter () {
		setIconImage (Starter.icon);
		Container root = getContentPane ();
		
		JPanel panelCodon = new JPanel (new BorderLayout ());
		JPanel panelEmbl = new JPanel (new BorderLayout ());
		
		radioCodon.setPreferredSize (new Dimension (120, 25));
		radioFasta.setPreferredSize (new Dimension (120, 25));
		radioEmbl.setPreferredSize (new Dimension (120, 25));
		
		fCodonFile.setPreferredSize (new Dimension (300, 25));
		fEmblFile.setPreferredSize (new Dimension (300, 25));
		fFastaFile.setPreferredSize (new Dimension (300, 25));
		fWorkspace.setPreferredSize (new Dimension (300, 25));
		
		panelFileChooseCodon.add (fCodonFile);
		
		panelFileChooseCodon.add (bOpenChooserCodon);
		panelCodon.add (radioCodon, BorderLayout.WEST);
		panelCodon.add (panelFileChooseCodon, BorderLayout.CENTER);
		
		panelFileChooseEmbl.add (fEmblFile);
		
		panelFileChooseEmbl.add (bOpenChooserEmbl);
		panelEmbl.add (radioEmbl, BorderLayout.WEST);
		panelEmbl.add (panelFileChooseEmbl, BorderLayout.CENTER);
		
		JPanel panelFasta = new JPanel (new BorderLayout ());
		panelFileChooseFasta.add (fFastaFile);
		panelFileChooseFasta.add (bOpenChooserFasta);
		JPanel panelOptions = new JPanel (new GridLayout (2, 1));
		panelOptions.add (ckStartBlastingQueue);
		ckStartBlastingQueue.setSelected (false);
		panelOptions.add (ckLoadEveryORFDetails);
		ckLoadEveryORFDetails.setSelected (true);
		panelFasta.add (radioFasta, BorderLayout.WEST);
		panelFasta.add (panelFileChooseFasta, BorderLayout.CENTER);
//		panelFasta.add (panelOptions, BorderLayout.SOUTH);
		
		ButtonGroup gp = new ButtonGroup ();
		gp.add (radioCodon);
		gp.add (radioFasta);
		gp.add (radioEmbl);
		radioFasta.setSelected (true);
			
		JPanel pint = new JPanel (new BorderLayout ());
		
		JPanel panelW = new JPanel (new BorderLayout ());
		JLabel l = new JLabel ("Workspace dir: ");
		l.setPreferredSize (new Dimension (115, 25));
		panelW.add (l, BorderLayout.WEST);
		panelW.add (fWorkspace);
		panelW.add (bOpenWorkspace, BorderLayout.EAST);
		panelW.add (new JPanel (), BorderLayout.SOUTH);
		
		pint.add (panelCodon, BorderLayout.NORTH);
		pint.add (panelEmbl, BorderLayout.CENTER);
		pint.add (panelFasta, BorderLayout.SOUTH);
		root.add (panelW, BorderLayout.NORTH);
		root.add (pint, BorderLayout.CENTER);
		root.add (bStart, BorderLayout.SOUTH);
		
		bOpenChooserFasta.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				FileNameExtensionFilter filter = new FileNameExtensionFilter ("FASTA", "fasta", "fa");
				chooser.setFileFilter (filter);
				chooser.setFileSelectionMode (JFileChooser.FILES_ONLY);
			    int returnVal = chooser.showOpenDialog (Starter.this);
			    if (returnVal == JFileChooser.APPROVE_OPTION) {
			    	try {
						fFastaFile.setText (chooser.getSelectedFile ().getCanonicalPath ());
					} catch (IOException e) {
						e.printStackTrace ();
					}
			    }
			}
		});
		bOpenChooserCodon.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				FileNameExtensionFilter filter = new FileNameExtensionFilter ("CODON", "codon");
				chooser.setFileFilter (filter);
				chooser.setFileSelectionMode (JFileChooser.FILES_ONLY);
			    int returnVal = chooser.showOpenDialog (Starter.this);
			    if (returnVal == JFileChooser.APPROVE_OPTION) {
			    	try {
						fCodonFile.setText (chooser.getSelectedFile ().getCanonicalPath ());
					} catch (IOException e) {
						e.printStackTrace ();
					}
			    }
			}
		});
		bOpenChooserEmbl.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				FileNameExtensionFilter filter = new FileNameExtensionFilter ("Embl & GenBank", "embl", "gb");
				chooser.setFileSelectionMode (JFileChooser.FILES_ONLY);
				chooser.setFileFilter (filter);
			    int returnVal = chooser.showOpenDialog (Starter.this);
			    if (returnVal == JFileChooser.APPROVE_OPTION) {
			    	try {
						fEmblFile.setText (chooser.getSelectedFile ().getCanonicalPath ());
					} catch (IOException e) {
						e.printStackTrace ();
					}
			    }
			}
		});
		bOpenWorkspace.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				chooser.setFileSelectionMode (JFileChooser.DIRECTORIES_ONLY);
			    int returnVal = chooser.showOpenDialog (Starter.this);
			    if (returnVal == JFileChooser.APPROVE_OPTION) {
			    	try {
						fWorkspace.setText (chooser.getSelectedFile ().getCanonicalPath ());
						chooser.setCurrentDirectory (chooser.getSelectedFile ());
					} catch (IOException e) {
						e.printStackTrace ();
					}
			    }
			}
		});
		radioCodon.addActionListener (new ActionListener () {
			public void actionPerformed(ActionEvent arg0) {
				setActive (1);
			}
		});
		radioFasta.addActionListener (new ActionListener () {
			public void actionPerformed(ActionEvent arg0) {
				setActive (0);
			}
		});
		radioEmbl.addActionListener (new ActionListener () {
			public void actionPerformed(ActionEvent arg0) {
				setActive (2);
			}
		});
		bStart.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent arg0) {
				String wDir = fWorkspace.getText ();
				File fwDir = new File (wDir);
				fwDir.mkdirs ();
				if (!fwDir.exists () || !fwDir.isDirectory ()) {
					JOptionPane.showMessageDialog (Starter.this, "Workspace must be e valid directory");
					return;
				}
				CodonConfiguration.workspace = wDir + "/";
				if (radioFasta.isSelected ()) {
					OrfRecognizer reco = getReconhecedor ();
					reco.load (fFastaFile.getText ());
					setVisible (false);
					saveConfig ();
					new PlotOrfsFrame (reco, ckStartBlastingQueue.isSelected (), ckLoadEveryORFDetails.isSelected (), ckLoadEveryORFDetails.isSelected ());
				}
				else if (radioCodon.isSelected ()) {
					OrfRecognizer reco = getReconhecedor ();
					File f = new File (fCodonFile.getText ());
					reco.loadFromCodonProject (f);
					setVisible (false);
					saveConfig ();
					PlotOrfsFrame p = new PlotOrfsFrame (reco, false, false, false);
					p.currentSaveFile = f;
					p.changed = false;
					p.updateFilter ();
				}
				else {
					OrfRecognizer reco = getReconhecedor ();
					reco.loadFromEmbl (new File (fEmblFile.getText ()));
					setVisible (false);
					saveConfig ();
					PlotOrfsFrame p = new PlotOrfsFrame (reco, false, true, false);
					p.updateFilter ();
					p.sendRNAS ();
				}
			}
		});
		setDefaultCloseOperation (DISPOSE_ON_CLOSE);
		setLocationByPlatform (true);
		setTitle ("CODON - Select the workspace directory and your entry file");
		pack ();
		setActive (0);
		loadConfig ();
		setVisible (true);
		setResizable (false);
		SwingUtilities.invokeLater (new Runnable () {
			public void run  () {
				bStart.grabFocus ();
			}
		});
	}
	
	public void setActive (int active) {
		fFastaFile.setEnabled (active == 0);
		bOpenChooserFasta.setEnabled (active == 0);
		ckLoadEveryORFDetails.setEnabled (active == 0);
		ckStartBlastingQueue.setEnabled (active == 0);
		fCodonFile.setEnabled (active == 1);
		bOpenChooserCodon.setEnabled (active == 1);
		fEmblFile.setEnabled (active == 2);
		bOpenChooserEmbl.setEnabled (active == 2);
	}
	
	public void loadConfig () {
		File f = new File ("data/config/init.cfg");
		if (!f.exists ()) return;
		try {
			BufferedReader reader = new BufferedReader (new FileReader (f));
			String l;
			while ((l = reader.readLine ()) != null) {
				String s [] = l.split ("=");
				if (s [0].equals ("project") && s.length > 1) fCodonFile.setText (s [1]);
				if (s [0].equals ("embl") && s.length > 1) fEmblFile.setText (s [1]);
				if (s [0].equals ("file") && s.length > 1) fFastaFile.setText (s [1]);
				if (s [0].equals ("workspace") && s.length > 1) {
					fWorkspace.setText (s [1]);
					chooser.setCurrentDirectory (new File (s [1]));
				}
				if (s [0].equals ("load_from") && s.length > 1) {
					int type = Integer.parseInt (s [1]);
					radioFasta.setSelected (type == 0);
					radioCodon.setSelected (type == 1);
					radioEmbl.setSelected (type == 2);
					setActive (type);
				}
				if (s [0].equals ("start_blasting") && s.length > 1) ckStartBlastingQueue.setSelected (Boolean.parseBoolean (s [1]));
				if (s [0].equals ("load_orfs") && s.length > 1) ckLoadEveryORFDetails.setSelected (Boolean.parseBoolean (s [1]));
			}
			reader.close ();
		} catch (IOException e) {
			e.printStackTrace ();
		}
	}
	
	public void saveConfig () {
		try {
			BufferedWriter writer = new BufferedWriter (new FileWriter ("data/config/init.cfg"));
			writer.write ("workspace=" + fWorkspace.getText ());
			writer.newLine ();
			writer.write ("project=" + fCodonFile.getText ());
			writer.newLine ();
			writer.write ("embl=" + fEmblFile.getText ());
			writer.newLine ();
			writer.write ("file=" + fFastaFile.getText ());
			writer.newLine ();
			writer.write ("load_from=" + (radioFasta.isSelected () ? 0 : (radioCodon.isSelected () ? 1 : 2)));
			writer.newLine ();
			writer.write ("start_blasting=" + ckStartBlastingQueue.isSelected ());
			writer.newLine ();
			writer.write ("load_orfs=" + ckLoadEveryORFDetails.isSelected ());
			writer.newLine ();
			writer.close ();
		} catch (IOException e) {
			e.printStackTrace ();
		}
	}
	
	public static OrfRecognizer getReconhecedor () {
		OrfRecognizer reco = new OrfRecognizer ();
		return reco;
	}
	
	public static void main (String [] args) {
//		try {
//			System.setOut (new PrintStream (new FileOutputStream ("log.txt")));
//		} catch (FileNotFoundException e) {
//			e.printStackTrace();
//		}
		try {
			// Set cross-platform Java L&F (also called "Metal")
//			UIManager.setLookAndFeel (UIManager.getCrossPlatformLookAndFeelClassName ());
//			UIManager.setLookAndFeel ("com.sun.java.swing.plaf.nimbus.NimbusLookAndFeel");
//			UIManager.setLookAndFeel ("com.sun.java.swing.plaf.windows.WindowsLookAndFeel");
			UIManager.setLookAndFeel ("javax.swing.plaf.nimbus.NimbusLookAndFeel");
		} 
		catch (UnsupportedLookAndFeelException e) {
			e.printStackTrace ();
		}
		catch (ClassNotFoundException e) {
			e.printStackTrace ();
		}
		catch (InstantiationException e) {
			e.printStackTrace ();
		}
		catch (IllegalAccessException e) {
			e.printStackTrace ();
		}
		CodonConfiguration.out = System.out;
		System.setOut (new PrintStream (new OutputStream () {
			public void write (int arg0) throws IOException {}
		}));
		new Starter ();
	}
	
}


