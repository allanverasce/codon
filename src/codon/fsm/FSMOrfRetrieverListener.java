package codon.fsm;

import java.util.ArrayList;

import codon.data.ORF;

public interface FSMOrfRetrieverListener {
	public void found (ArrayList <ORF> orfFound); 
}
