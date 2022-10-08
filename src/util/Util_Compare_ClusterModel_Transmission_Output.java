package util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.compress.archivers.sevenz.SevenZArchiveEntry;
import org.apache.commons.compress.archivers.sevenz.SevenZFile;

import sim.Simulation_ClusterModelTransmission;

public class Util_Compare_ClusterModel_Transmission_Output {

	File ref_result_dir;
	File results_dir;
	int[][] analysis_col_select;

	final static int BUFFER = 2048;
	final static String USAGE_INFO = String.format(
			"Usage: java %s -compare RESULTS_DIRECTORY REF_RESULT_DIRNAME ANALSIS_COL_SELECT",
			Util_Compare_ClusterModel_Transmission_Output.class.getName());

	final static int[] ANALYSIS_TYPE_OPTIONS = new int[] {
			Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE };

	public Util_Compare_ClusterModel_Transmission_Output(File results_dir, File ref_result_dir,
			int[][] analysis_col_select) {
		super();
		this.results_dir = results_dir;
		this.ref_result_dir = ref_result_dir;
		this.analysis_col_select = analysis_col_select;
	}

	public void startCompare() throws IOException {
		File[] compare_result_dir = results_dir.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.isDirectory() && !pathname.getName().equals(ref_result_dir.getName());
			}
		});

		Arrays.sort(compare_result_dir, new Comparator<File>() {
			@Override
			public int compare(File o1, File o2) {
				return o1.getName().compareTo(o2.getName());
			}

		});

		System.out.printf("Comparing %d result(s) from %s, using %s as reference.\n", compare_result_dir.length,
				results_dir.getAbsolutePath(), ref_result_dir.getName());

		Pattern patten_src_zip_file;
		int[] col_sel;

		for (int analysis_type : ANALYSIS_TYPE_OPTIONS) {
			col_sel = new int[0]; // Length of zero = all, null = non selected

			if (analysis_type < analysis_col_select.length) {
				col_sel = analysis_col_select[analysis_type];
			}

			if (col_sel != null) {
				switch (analysis_type) {
				case Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE:
					// Prevalence by person
					patten_src_zip_file = Pattern
							.compile(Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_PERSON_ZIP.replaceAll("%d",
									"(-{0,1}(?!0)\\\\d+)"));

					File[] src_zips = ref_result_dir.listFiles(generateFileFilterByPattern(patten_src_zip_file));

					HashMap<String, ArrayList<String[]>> zip_ent = new HashMap<>();

					for (File zipFile : src_zips) {
						Matcher m;
						m = patten_src_zip_file.matcher(zipFile.getName());
						m.matches();
						long cMapSeed = Long.parseLong(m.group(1));
						zip_ent = extractedLinesFrom7Zip(zipFile, zip_ent);
					}
					
					
					// TODO: Match with other mapping in compare_result_dir

					break;
				default:
					System.out.printf("Analysis type %d not support yet\n", analysis_type);

				}
			}

		}

	}

	private static HashMap<String, ArrayList<String[]>> extractedLinesFrom7Zip(File zipFile,
			HashMap<String, ArrayList<String[]>> zip_ent) throws IOException {
		SevenZFile inputZip = new SevenZFile(zipFile);
		SevenZArchiveEntry inputEnt;

		byte[] buf = new byte[BUFFER];
		while ((inputEnt = inputZip.getNextEntry()) != null) {
			String file_name = inputEnt.getName();
			StringBuilder str_builder = new StringBuilder();
			String line;
			ArrayList<String[]> lines = new ArrayList<>();
			int count;
			while ((count = inputZip.read(buf, 0, BUFFER)) != -1) {
				str_builder.append(new String(Arrays.copyOf(buf, count)));
			}
			BufferedReader reader = new BufferedReader(new StringReader(str_builder.toString()));
			while ((line = reader.readLine()) != null) {
				lines.add(line.split(","));
			}
			zip_ent.put(file_name, lines);
		}
		inputZip.close();
		return zip_ent;
	}

	private static FileFilter generateFileFilterByPattern(Pattern pattern) {
		FileFilter filter = new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pattern.matcher(pathname.getName()).matches();
			}
		};

		return filter;
	}

	public static void launch(String[] arg) throws IOException {
		if (arg.length < 3) {
			System.out.println(Util_Compare_ClusterModel_Transmission_Output.USAGE_INFO);
			System.exit(1);
		} else {
			File results_dir = new File(arg[0]);
			File ref_result_dir = new File(results_dir, arg[1]);
			int[][] analysis_col_select = (int[][]) PropValUtils.propStrToObject(arg[2], int[][].class);

			Util_Compare_ClusterModel_Transmission_Output cmp = new Util_Compare_ClusterModel_Transmission_Output(
					results_dir, ref_result_dir, analysis_col_select);

			cmp.startCompare();
		}

	}

}
