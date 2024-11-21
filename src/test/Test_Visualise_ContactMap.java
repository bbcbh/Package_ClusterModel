package test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;

import relationship.ContactMap;
import sim.Simulation_ClusterModelGeneration;
import visualise.Visualise_ContactMap;

public class Test_Visualise_ContactMap {

	public static void main(String[] args) throws IOException {

		final long cMap_seed = -1313307668265488683l;
		final int time_max = 10950;
		final int[] popComposition = new int[] { 0, 0, 14000, 0 };		
		final File baseDir = new File("C:\\Users\\bhui\\OneDrive - UNSW\\Bridging_model\\GenMap_BALI");
		
		long tic;
		
		// Reading of cMap		
		File cMapFile = new File(baseDir,
				String.format(Simulation_ClusterModelGeneration.FILENAME_FORMAT_ALL_CMAP, cMap_seed, time_max));
		
		
		System.out.printf("Reading cMap CSV from %s ...\n", cMapFile.getAbsolutePath());
		
		tic = System.currentTimeMillis();				
		StringWriter cMap_str = new StringWriter();
		PrintWriter pWri = new PrintWriter(cMap_str);
		
		BufferedReader reader = new BufferedReader(new FileReader(cMapFile));
		String line;		
		while ((line = reader.readLine()) != null) {
			pWri.println(line);
		}
		pWri.close();
		reader.close();		
		
		System.out.printf("Reading CSV completed. Time reqired = %.3fs\n" , (System.currentTimeMillis() - tic)/ 1000f);
		
		tic = System.currentTimeMillis();
		ContactMap cMap = ContactMap.ContactMapFromFullString(cMap_str.toString());
		System.out.printf("Generation of cMap completed. Time reqired = %.3fs\n" , (System.currentTimeMillis() - tic)/ 1000f);
		
		// Visualise
		tic = System.currentTimeMillis();
		Visualise_ContactMap.visualiseCluster(baseDir, popComposition, cMap, cMap_seed, 1,
				Visualise_ContactMap.FLAG_GEN_PNG);		
		System.out.printf("Visualisation of cMap completed. Time reqired = %.3fs\n" , (System.currentTimeMillis() - tic)/ 1000f);

	}

}
