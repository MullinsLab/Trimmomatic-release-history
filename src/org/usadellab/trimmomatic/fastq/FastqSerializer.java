package org.usadellab.trimmomatic.fastq;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.zip.GZIPOutputStream;

public class FastqSerializer {

	private PrintStream stream;
	
	public FastqSerializer()
	{
		
	}
	
	public void open(File file) throws IOException
	{
		String name=file.getName();
	
		if(name.endsWith(".gz"))
			stream=new PrintStream(new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(file)),1000000));
		else
			stream=new PrintStream(new BufferedOutputStream(new FileOutputStream(file),1000000));
	}

	public void open(PrintStream stream) throws IOException
	{
		this.stream=new PrintStream(stream);
	}
	
	public void close() throws IOException
	{
		stream.close();
	}
	
	public synchronized void writeRecord(FastqRecord record) 
	{
		stream.println("@"+record.getName());
		stream.println(record.getSequence());
		stream.println("+"+record.getComment());
		stream.println(record.getQuality());
	}
	
}
