package test;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import util.Util_7Z_CSV_Entry_Extract_Callable;
import util.Util_SimulationDirectoryModifications;

public class Test_PreSim_GenerateDirWithPredefineSeed_Bali {
	public static void main(String[] args) throws IOException {
		File seedDir = new File("C:\\Users\\bhui\\OneDrive - UNSW\\Results_Bali\\"
				+ "20240824 (Poster)\\Bali_Baseline\\SimClusterModel_Transmission");
		File destDir = new File("C:\\Users\\bhui\\Documents\\Java_Test\\Result_Bali\\Bali_Baseline");
		File sourcePBS = new File(
				"C:\\Users\\Bhui\\Documents\\Java_Test\\Prop_Template\\PBS\\ClusterModel_Sim_Template_Bali.pbs");

		String folderFormat = "SimClusterModel_Transmission_%03d";
		String seedFileName = "Seed_List.csv";
		int numSimPerDir = 8;
		
		String pbsIdFormat = "#PBS -N BaliU100_%03d";
		String pbs_filename_format = "ClusterModel_Sim_%03d.pbs";

		File[] directCopyFiles = new File[] { 
				new File(destDir, "simSpecificSim.prop"),
				new File(destDir, "simSpecificSwitch.prop"), 
				new File(destDir, "dx_Bali.prop"), 
				};

		File[] seedFiles = new File[] { new File(seedDir, "Seed_List.csv"), new File(seedDir, "Seed_List.csv.zip"), };

		// Extract seed file
		ArrayList<String[]> seedFileTxtArray = Util_SimulationDirectoryModifications.extractSeedFilesEntries(seedFiles);

		System.out.printf("%d seed file extracted from %s.\n", seedFileTxtArray.size(), seedDir.getAbsolutePath());

		Collections.sort(seedFileTxtArray, new Comparator<String[]>() {
			@Override
			public int compare(String[] o1, String[] o2) {
				return o1[1].compareTo(o2[1]);
			}
		});

		// Creating folder and associated PBS
		int folderNum = 0;
		
		File pbs_folder = new File(destDir.getParent(), "PBS_" + destDir.getName());
		pbs_folder.mkdirs();

		PrintWriter batchScriptWriter = new PrintWriter(new File(pbs_folder, "batch_qsub")) {
			@Override
			public void println() {
				write('\n');
			}
		};

		batchScriptWriter.println("#!/bin/bash");
		batchScriptWriter.println("echo \"Batch submiting...\"");

		String[] line_pbs_arr = Util_7Z_CSV_Entry_Extract_Callable.extracted_lines_from_text(sourcePBS);

		for (String[] srcTxt : seedFileTxtArray) {
			File simDir = new File(destDir, String.format(folderFormat, folderNum));
			simDir.mkdirs();

			// Seed file
			PrintWriter seedWri = new PrintWriter(new File(simDir, seedFileName));
			for (String str : srcTxt) {
				seedWri.println(str);
			}
			seedWri.close();

			// Direct copy files
			for (File directCopy : directCopyFiles) {
				File newDest = new File(simDir, directCopy.getName());
				Files.copy(directCopy.toPath(), newDest.toPath());
			}

			// PBS file
			File pbsFile = new File(pbs_folder, String.format(pbs_filename_format, folderNum));

			PrintWriter pbs_Writer = new PrintWriter(pbsFile) {
				@Override
				public void println() {
					write('\n');
				}
			};

			for (int n = 0; n < line_pbs_arr.length; n++) {
				String lineEnt = line_pbs_arr[n];
				switch (n) {
				case 5:
					lineEnt = String.format(pbsIdFormat, folderNum);
					break;
				case 9:
					lineEnt = String.format("#PBS -l ncpus=%d", Math.min(numSimPerDir, srcTxt.length - 1));
					break;
				case 14:
					lineEnt = String.format("#PBS -o ./%s.out", pbsFile.getName());
					break;
				case 15:
					lineEnt = String.format("#PBS -e ./%s.err", pbsFile.getName());
					break;
				case 23:
					lineEnt = String.format(
							"java -jar \"ClusterModel.jar\" "
									+ "-trans \"/srv/scratch/z2251912/All_Transmission/%s/%s\" "
									+ "-export_skip_backup -printProgress -seedMap=Seed_List.csv",
							destDir.getName(), simDir.getName());
				default:

				}
				pbs_Writer.println(lineEnt);
			}
			pbs_Writer.close();

			batchScriptWriter.printf("qsub %s\n", pbsFile.getName());

			folderNum++;
		}
		
		batchScriptWriter.println("echo \"Submit completed!\"");
		batchScriptWriter.close();

		System.out.printf("%d seed file directories generated at %s.\n", folderNum, destDir.getAbsolutePath());

	}

}
