package org.usadellab.trimmomatic.fastq;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;

public class FastqParser
{
	private int phredOffset;
	
	private BufferedReader reader;

	private FastqRecord current;

	public FastqParser(int phredOffset)
	{
		this.phredOffset=phredOffset;
	}

	public void parseOne() throws IOException
	{
		current = null;

		String name;
		String sequence;
		String comment;
		String quality;
		
		String line;
		
		line=reader.readLine();
		if(line==null)
			return;
		
		if(line.startsWith("@"))
			name=line.substring(1);
		else
			throw new RuntimeException("Invalid FASTQ name line: "+line);
		
		sequence=reader.readLine();
	
		line=reader.readLine();
		if(line.startsWith("+"))
			comment=line.substring(1);
		else
			throw new RuntimeException("Invalid FASTQ comment line: "+line);

		quality=reader.readLine();

		current=new FastqRecord(name, sequence, comment, quality, phredOffset);
	}

	public void parse(File file) throws IOException
	{
		String name=file.getName();
		
		if(name.endsWith(".gz"))
			reader=new BufferedReader(new InputStreamReader(new GZIPInputStream(new BufferedInputStream(new FileInputStream(file),1000000))));
		else			
			reader = new BufferedReader(new InputStreamReader(new BufferedInputStream(new FileInputStream(file),1000000)));
		
		parseOne();
	}
	
	public void close() throws IOException
	{
		reader.close();
	}

	public boolean hasNext()
	{
		return current != null;
	}

	public FastqRecord next() throws IOException
	{
		FastqRecord current = this.current;
		parseOne();

		return current;
	}
}
