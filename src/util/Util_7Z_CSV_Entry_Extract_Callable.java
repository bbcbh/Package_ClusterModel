package util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.StringReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.compress.archivers.sevenz.SevenZArchiveEntry;
import org.apache.commons.compress.archivers.sevenz.SevenZFile;

public class Util_7Z_CSV_Entry_Extract_Callable implements Callable<Map<String, float[]>> {

	final private File BaseDir;
	final private String[] CSV_Filenames_Pattern;
	final private String CMapId;
	final private int[][] Entry_Pairs; // {[File_id, colNum}

	public Util_7Z_CSV_Entry_Extract_Callable(File baseDir, String cMapId, String[] csv_Filenames_Pattern, int[][] entry_Pairs) {
		super();
		BaseDir = baseDir;
		CMapId = cMapId;
		CSV_Filenames_Pattern = csv_Filenames_Pattern;
		Entry_Pairs = entry_Pairs;
	}

	@Override
	public Map<String, float[]> call() throws Exception {
		final int BUFFER = 2048;
		final byte[] buf = new byte[BUFFER];
		HashMap<String, float[]> res = new HashMap<>();

		for (int file_id = 0; file_id < CSV_Filenames_Pattern.length; file_id++) {

			Pattern pattern_archive_file = Pattern
					.compile((CSV_Filenames_Pattern[file_id].replaceFirst("_%d", "") + ".7z").replaceAll("%d",CMapId));
			Pattern pattern_csv_file = Pattern
					.compile(CSV_Filenames_Pattern[file_id].replaceAll("%d", "(-{0,1}(?!0)\\\\d+)"));

			File[] archiveFile_Arr = BaseDir.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return pattern_archive_file.matcher(pathname.getName()).matches();
				}
			});

			for (File archiveFile : archiveFile_Arr) {
				SevenZFile archive7Z = new SevenZFile(archiveFile);
				SevenZArchiveEntry ent;
				while ((ent = archive7Z.getNextEntry()) != null) {
					Matcher m = pattern_csv_file.matcher(ent.getName());
					if (m.matches()) {
						String res_key = String.format("%s_%s", m.group(1), m.group(2));
						float[] res_ent = res.get(res_key);
						if (res_ent == null) {
							res_ent = new float[Entry_Pairs.length];
							res.put(res_key, res_ent);
						}
						StringBuilder txt_entries = new StringBuilder();
						int count;
						while ((count = archive7Z.read(buf, 0, BUFFER)) != -1) {
							txt_entries.append(new String(Arrays.copyOf(buf, count)));
						}
						String line, lastline = "";
						BufferedReader lines = new BufferedReader(new StringReader(txt_entries.toString()));
						while ((line = lines.readLine()) != null) {
							if (line.length() > 0) {
								lastline = line;
							}
						}
						lines.close();
						String[] lastlineEnt = lastline.split(",");
						for (int k = 0; k < Entry_Pairs.length; k++) {
							if (Entry_Pairs[k][0] == file_id) {
								res_ent[k] = Float.parseFloat(lastlineEnt[Entry_Pairs[k][1]]);
							}
						}
					}
				}
				archive7Z.close();
			}
		}
		return res;
	}

}
