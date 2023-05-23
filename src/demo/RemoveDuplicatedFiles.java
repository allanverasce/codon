package demo;

import java.io.File;

public class RemoveDuplicatedFiles {
	public static void main(String[] args) {
		File f = new File ("data/ERR2348864/blasted");
		File lst [] = f.listFiles ();
		for (int i = 0; i < lst.length; i++) {
			if (lst [i].getName ().indexOf ("(") != -1) lst [i].delete ();
			if (lst [i].length() == 0) lst [i].delete ();
		}
	}
}
