package codon.ui;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.JToggleButton;
import javax.swing.SpinnerNumberModel;

import codon.data.ORF;
import codon.recognizer.OrfRecognizer;

public class FilterBox extends JDialog {
	private static final long serialVersionUID = 1L;
	
	static public class UserFilter {
		String name = "";
		ArrayList <VirtualFilterCel> cells = new ArrayList <> ();
		int action = 0;
		public UserFilter (String name, ArrayList <VirtualFilterCel> cells, int action) {
			this.name = name;
			this.cells = cells;
			this.action = action;
		}
	}
	
	static ArrayList <UserFilter> savedFilters = new ArrayList <> ();
	
	static class VirtualFilterCel {
		boolean not;
		int operation;
		int value;
		VirtualFilterCel (boolean not, int operation, int value) {
			this.not = not;
			this.operation = operation;
			this.value = value;
		}
		public boolean check (ORF orf) {
			if (operation == 0) {
				if (orf.getAccuracy () < value) return !not;
				else return not;
			}
			else if (operation == 1) {
				if (orf.getIdentity () < value) return !not;
				else return not;
			}
			else if (operation == 2) {
				if (orf.length < value) return !not;
				else return not;
			}
			else if (operation == 3) {
				if (orf.isGene ()) return !not;
				else return not;
			}
			else if (operation == 4) {
				if (!orf.isGene () && !orf.isUncharacterized ()) return !not;
				else return not;
			}
			else if (operation == 5) {
				if (orf.isUncharacterized ()) return !not;
				else return not;
			}
			else if (operation == 6) {
				if (!orf.isBlasted ()) return !not;
				else return not;
			}
			return false;
		}
	}

	static public void loadFilters () {
		File fFile = new File ("data/config/filters.txt");
		if (!fFile.exists ()) return;
		try {
			BufferedReader reader = new BufferedReader (new FileReader (fFile));
			String line;
			while ((line = reader.readLine ()) != null) {
				String name = line;
				int action = Integer.parseInt (reader.readLine ());
				ArrayList <VirtualFilterCel> cells = new ArrayList <> ();
				while (!(line = reader.readLine ()).equals ("<")) {
					String fields [] = line.split (";");
					VirtualFilterCel cel = new VirtualFilterCel (Boolean.parseBoolean (fields [0]),	Integer.parseInt (fields [1]), Integer.parseInt (fields [2]));
					cells.add (cel);
				}
				savedFilters.add (new UserFilter (name, cells, action));
			}
			
			reader.close ();
		} catch (IOException e) {
			e.printStackTrace ();
		}
	}
	
	static {
		loadFilters ();
	}
	
	static public void saveFilters () {
		try {
			BufferedWriter writer = new BufferedWriter (new FileWriter ("data/config/filters.txt"));
			for (UserFilter filter: savedFilters) {
				writer.write (filter.name); 
				writer.newLine ();
				writer.write (filter.action + "");
				writer.newLine ();
				for (VirtualFilterCel cel: filter.cells) {
					writer.write (cel.not + ";" + cel.operation + ";" + cel.value);
					writer.newLine ();
				}
				writer.write ("<");
				writer.newLine ();
			}
			writer.close ();
		} catch (IOException e) {
			e.printStackTrace ();
		}
	}
	
	static void addFilter (String name, ArrayList <FilterCel> cells, int operation) {
		for (UserFilter uF: savedFilters) {
			if (uF.name.equals (name)) {
				savedFilters.remove (uF);
				break;
			}
		}
		ArrayList <VirtualFilterCel> vcells = new ArrayList <FilterBox.VirtualFilterCel> ();
		for (FilterCel fC : cells) {
			boolean not = fC.getNot ();
			String filter = fC.getFilter ();
			int fOp = 0;
			int value = fC.getValue ();
			for (int i = 0; i < items.length; i++) {
				if (items[i].equals (filter)) {
					fOp = i;
					break;
				}
			}
			vcells.add (new VirtualFilterCel (not, fOp, value));
		}
		UserFilter uF = new UserFilter (name, vcells, operation);
		savedFilters.add (uF);
		saveFilters ();
	}
	
	static void applyFilter (PlotOrfsFrame frame, String name) {
		UserFilter filter = null;
		for (UserFilter uF: savedFilters) {
			if (uF.name.equals (name)) {
				filter = uF;
				break;
			}
		}
		if (filter == null) {
			JOptionPane.showMessageDialog (frame, "ERROR: Unknown filter " + name);
			return;
		}
		boolean remove = false;
		boolean add = false;
		
		if (filter.action == 1) add = true;
		else if (filter.action == 2) remove = true;
		else if (filter.action == 3) {
			remove = true;
			add = true;
		}
		
		for (ArrayList <ORF> orfs : frame.reco.orfsPorFrame) {
			for (ORF orf: orfs) {
				boolean in = true;
				if (add && orf.isRemoved () == remove) continue;
				for (int i = 0; i < filter.cells.size (); i++) {
					in &= filter.cells.get (i).check (orf);
				}
				if (in) orf.setRemoved (remove);
				else orf.setRemoved (!remove);
			}
		}
	}
		
	static String [] items = {
		"accuracy < ",
		"identity < ",
		"orf length < ",
		"gene",
		"other product",
		"uncharacterized ",
		"unblasted ",
	};
	
	class FilterCel extends JPanel {
		private static final long serialVersionUID = 1L;
		
		JLabel lInclude = new JLabel ("AND ");
		JToggleButton toggleNot = new JToggleButton (" has ");
		JComboBox <String> cSelect = new JComboBox <> (items);
		JSpinner spVal = new JSpinner ();
		JButton removeButton = new JButton ("x");
		SpinnerNumberModel smPercent = new SpinnerNumberModel (80, 0, 100, 1);
		SpinnerNumberModel smLength = new SpinnerNumberModel (80, 0, 100000, 10);
		
		public void setFirst (boolean first) {
			lInclude.setText (first ? "" : "AND ");
		}
		
		public FilterCel () {
			add (lInclude);
			add (toggleNot);
			add (cSelect);
			add (spVal);
			add (removeButton);
			spVal.setModel (smPercent);
			lInclude.setPreferredSize(new Dimension (40, 20));
			toggleNot.setPreferredSize(new Dimension (70, 20));
			spVal.setPreferredSize(new Dimension (90, 20));
			toggleNot.setMargin (new Insets (2, 2, 2, 2));
			toggleNot.addActionListener (new ActionListener () {
				public void actionPerformed (ActionEvent e) {
					if (toggleNot.isSelected ()) {
						if (toggleNot.getText ().equals (" has ")) {
							toggleNot.setText (" has not ");
						}
						else {
							toggleNot.setText (" is not ");
						}
					}
					else {
						if (toggleNot.getText ().equals (" has not ")) {
							toggleNot.setText (" has ");
						}
						else {
							toggleNot.setText (" is ");
						}
					}
				}
			});
			cSelect.addItemListener (new ItemListener () {
				public void itemStateChanged (ItemEvent e) {
					if (e.getStateChange () == ItemEvent.SELECTED) {
						String sel = (String) e.getItem ();
						if (sel == items [0]) {
							spVal.setEnabled (true);
							toggleNot.setText (toggleNot.isSelected () ? " has not " : " has ");
							spVal.setModel (smPercent);
						}
						else if (sel == items [1]) {
							spVal.setEnabled (true);
							toggleNot.setText (toggleNot.isSelected () ? " has not " : " has ");
							spVal.setModel (smPercent);
						} 
						else if (sel == items [2]) {
							spVal.setEnabled (true);
							toggleNot.setText (toggleNot.isSelected () ? " is not " : " is ");
							spVal.setModel (smLength);
						}
						else if (sel == items [3]) {
							spVal.setEnabled (false);
							toggleNot.setText (toggleNot.isSelected () ? " is not " : " is ");
						}
						else if (sel == items [4]) {
							spVal.setEnabled (false);
							toggleNot.setText (toggleNot.isSelected () ? " is not " : " is ");
						}
						else if (sel == items [5]) {
							spVal.setEnabled (false);
							toggleNot.setText (toggleNot.isSelected () ? " is not " : " is ");
						}
						else if (sel == items [6]) {
							spVal.setEnabled (false);
							toggleNot.setText (toggleNot.isSelected () ? " is not " : " is ");
						}
					}
				}
			});
			removeButton.addActionListener (new ActionListener () {
				public void actionPerformed (ActionEvent e) {
					cells.remove (FilterCel.this);
					fillPanelFilters ();
					pack ();
				}
			});
		}
		public boolean check (ORF orf) {
			String sel = (String) cSelect.getSelectedItem ();
			int v = (Integer) spVal.getValue ();
			boolean not = toggleNot.isSelected ();
			if (sel == items [0]) {
				if (orf.getAccuracy () < v) return !not;
				else return not;
			}
			else if (sel == items [1]) {
				if (orf.getIdentity () < v) return !not;
				else return not;
			}
			else if (sel == items [2]) {
				if (orf.length < v) return !not;
				else return not;
			}
			else if (sel == items [3]) {
				if (orf.isGene ()) return !not;
				else return not;
			}
			else if (sel == items [4]) {
				if (!orf.isGene () && !orf.isUncharacterized ()) return !not;
				else return not;
			}
			else if (sel == items [5]) {
				if (orf.isUncharacterized ()) return !not;
				else return not;
			}
			else if (sel == items [6]) {
				if (!orf.isBlasted ()) return !not;
				else return not;
			}
			return false;
		}
		public String getFilter () {
			return (String) cSelect.getSelectedItem ();
		}
		public boolean getNot () {
			return toggleNot.isSelected ();
		}
		public int getValue () {
			return (Integer) spVal.getValue ();
		}
	}
	
	JPanel panelFilters = new JPanel ();
	JPanel panelAdd = new JPanel ();
	ArrayList <FilterCel> cells = new ArrayList <FilterBox.FilterCel> ();
	OrfRecognizer reco;
	PlotOrfsFrame frame;
	JTextField fSave = new JTextField ();
	JCheckBox ckSave = new JCheckBox ("Save filter as: ");
	
	public FilterBox (PlotOrfsFrame frame, OrfRecognizer reco) {
		setTitle ("Filter");
		this.reco = reco;
		this.frame = frame;
		Container root = getContentPane ();
		
		JPanel savePanel = new JPanel (new BorderLayout ());
		fSave.setEnabled (false);
		savePanel.add (ckSave, BorderLayout.WEST);
		savePanel.add (fSave);
		root.add (savePanel, BorderLayout.NORTH);
		ckSave.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent e) {
				fSave.setEnabled (ckSave.isSelected ());
			}
		});
		
		root.add (panelFilters);
		cells.add (new FilterCel ());
		fillPanelFilters ();
		JButton add = new JButton ("+");
		panelAdd.add (add);
		add.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent e) {
				cells.add (new FilterCel ());
				fillPanelFilters ();
				pack ();
			}
		});
		JPanel actionPane = new JPanel (new GridLayout (1, 4));
		JButton bAddV = new JButton ("Add Visible");
		JButton bSetV = new JButton ("Set Visible");
		JButton bAddR = new JButton ("Add Removed");
		JButton bSetR = new JButton ("Set Removed");
		actionPane.add (bAddV);
		actionPane.add (bSetV);
		actionPane.add (bAddR);
		actionPane.add (bSetR);
		bAddV.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent e) {
				applyFilter (false, true);
			}
		});
		bSetV.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent e) {
				applyFilter (false, false);
			}
		});
		bAddR.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent e) {
				applyFilter (true, true);
			}
		});
		bSetR.addActionListener (new ActionListener () {
			public void actionPerformed (ActionEvent e) {
				applyFilter (true, false);
			}
		});
		root.add (actionPane, BorderLayout.SOUTH);
		pack ();
		setLocationRelativeTo (frame);
		setDefaultCloseOperation (DISPOSE_ON_CLOSE);
		setVisible (true);
	}
	
	public void applyFilter (boolean remove, boolean add) {
		for (ArrayList <ORF> orfs : reco.orfsPorFrame) {
			for (ORF orf: orfs) {
				boolean in = true;
				if (add && orf.isRemoved () == remove) continue;
				for (int i = 0; i < cells.size (); i++) {
					in &= cells.get (i).check (orf);
				}
				if (in) orf.setRemoved (remove);
				else orf.setRemoved (!remove);
			}
		}
		frame.updateFilter ();
		if (ckSave.isSelected ()) {
			int operation = 0;
			if (!remove && add) operation = 1;
			if (remove && !add) operation = 2;
			if (remove && add) operation = 3;
			addFilter (fSave.getText (), cells, operation);
			frame.buildFilterList ();
		}
	}
	
	public void fillPanelFilters () {
		panelFilters.removeAll ();
		panelFilters.setLayout (new GridLayout (cells.size () + 1, 1));
		for (FilterCel f: cells) {
			panelFilters.add (f);
			f.setFirst (false);
		}
		panelFilters.add (panelAdd);
		if (cells.size () > 0) cells.get (0).setFirst (true);
	}
	
}
