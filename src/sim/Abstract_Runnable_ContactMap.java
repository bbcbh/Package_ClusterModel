package sim;

import java.io.File;

public abstract class Abstract_Runnable_ContactMap implements Runnable {
	
	protected String runnableId = null;
	protected File baseDir = null;
	
	
	public String getRunnableId() {
		return runnableId;
	}

	public void setRunnableId(String runnableId) {
		this.runnableId = runnableId;
	}

	public void setBaseDir(File baseDir) {
		this.baseDir = baseDir;
	}

	public abstract Object[] getRunnable_fields();
	

}
