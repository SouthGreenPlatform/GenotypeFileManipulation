/*******************************************************************************
 * Genotype file manipulation - Helper for reading PLINK and Eigenstrat files
 * Copyright (C) 2018, <CIRAD>
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License, version 3 as published by
 * the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
 * details.
 *
 * See <http://www.gnu.org/licenses/agpl.html> for details about GNU General
 * Public License V3.
 *******************************************************************************/
package fr.cirad.tools.genotypes;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.LineNumberReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.StringTokenizer;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import org.apache.commons.io.FileUtils;
import org.apache.log4j.Logger;

import com.sun.org.apache.xpath.internal.functions.WrongNumberArgsException;

public class PlinkEigenstratTool {
	
	private static final Logger LOG = Logger.getLogger(PlinkEigenstratTool.class);
		
	static public HashMap<String /*preferred Illumina ID*/, List<String> /*Illumina synonyms*/> readVariantSynonyms(File variantSynonymFile) throws FileNotFoundException
	{
		HashMap<String, List<String>> variantSynonyms = new HashMap<String, List<String>>();
		Scanner refVariantScanner = new Scanner(variantSynonymFile);
		while (refVariantScanner.hasNextLine())
		{
			String[] splittedLine = refVariantScanner.nextLine().split("=");
			if (splittedLine.length != 2)
				continue;
			
			variantSynonyms.put(splittedLine[0], Arrays.asList(splittedLine[1].split(",")));
		}
		refVariantScanner.close();
		return variantSynonyms;
	}
	
	static public HashMap<String /*variant id or synonym*/, HashMap<String /*PLink genotype code*/, Integer /*corresponding Eigenstrat genotype code*/>> readVariantGenotypeCodes(File genotypeCodeConversionFile) throws FileNotFoundException
	{
		HashMap<String /*variant id or synonym*/, HashMap<String /*PLink genotype code*/, Integer /*corresponding Eigenstrat genotype code*/>> perVariantGtCodes = new HashMap<String, HashMap<String, Integer>>();		
		Scanner scanner = new Scanner(genotypeCodeConversionFile);
		while (scanner.hasNextLine())
		{	// looping on ref variants
			String sLine = scanner.nextLine();
			String[] synonymAndData = sLine.split("\t");
			if (synonymAndData.length < 2)
				continue;
			
			HashMap<String, Integer> synonymDataMap = new HashMap<String, Integer>();
			for (int i=1; i<synonymAndData.length; i++)
			{
				String[] genotypeMapping = synonymAndData[i].split(":");
				if (genotypeMapping.length < 2)
					continue;

				synonymDataMap.put(genotypeMapping[0], Integer.parseInt(genotypeMapping[1]));
			}
			perVariantGtCodes.put(synonymAndData[0], synonymDataMap);
		}
		scanner.close();
		return perVariantGtCodes;
	}
	
	
	static public LinkedHashMap<String, String> getVariantsAndPositionsFromPlinkMapFile(File mapFile, LinkedHashSet<Integer> redundantVariantIndexesToFill, String sSequenceAndPositionSeparator) throws Exception
	{
		return getVariantsAndPositionsFromPlinkMapFile(mapFile.toURI().toURL(), redundantVariantIndexesToFill, sSequenceAndPositionSeparator);
	}

	static public LinkedHashMap<String, String> getVariantsAndPositionsFromPlinkMapFile(URL mapFileUrl, LinkedHashSet<Integer> redundantVariantIndexesToFill, String sSequenceAndPositionSeparator) throws Exception
	{
		StringBuffer errors = new StringBuffer();
		
		LinkedHashMap<String, String> variantsInMapFile = new LinkedHashMap<String, String>();
		Scanner userMapScanner = new Scanner(mapFileUrl.openStream());
		int nCurrentLine = -1;
		while (userMapScanner.hasNextLine())
		{
			nCurrentLine++;
			String sLine = userMapScanner.nextLine().replaceAll("\t", " ");
			if (sLine.trim().length() == 0)
			{
				errors.append("\n- Found empty line in .map file at position " + nCurrentLine);
				continue;
			}

			StringTokenizer st = new StringTokenizer(sLine, " ");
			String sSeq = st.nextToken();
			String sVariant = st.nextToken();
			st.nextToken(); // skip Genetic distance in morgans
			String sBpPos = st.nextToken();
			if (variantsInMapFile.containsKey(sVariant))
			{
				LOG.warn("Variant " + sVariant + " is defined several times in MAP file");
				redundantVariantIndexesToFill.add(nCurrentLine);
			}
			variantsInMapFile.put(sVariant, sSeq + sSequenceAndPositionSeparator + sBpPos);
		}
		
		userMapScanner.close();
		if (errors.length() > 0)
			throw new Exception("Uploaded data is invalid: \n" + errors.toString());

		return variantsInMapFile;
	}
	
	static public String[] getVariantsFromPlinkMapFile(File mapFile, LinkedHashSet<Integer> redundantVariantIndexesToFill) throws Exception
	{
		LinkedHashMap<String, String> result = getVariantsAndPositionsFromPlinkMapFile(mapFile, redundantVariantIndexesToFill, "::");
		return result.keySet().toArray(new String[result.size()]);
	}

	static public ArrayList<String> getVariantsFromEigenstratSnpFile(File snpFile) throws FileNotFoundException
	{
		ArrayList<String> refVariants = new ArrayList<String>();
		Scanner refSnpScanner = new Scanner(snpFile);
		while (refSnpScanner.hasNextLine())
		{
			String sLine = refSnpScanner.nextLine();
			refVariants.add(sLine.substring(0, sLine.indexOf("\t")));
		}
		refSnpScanner.close();
		return refVariants;
	}

	static public LinkedHashMap<String, String> getIndividualToPopulationMapFromEigenstratIndFile(File indFile) throws FileNotFoundException
	{
		LinkedHashMap<String /*ref individual*/, String /*population*/> refIndividualToPopulationMap = new LinkedHashMap<String, String>();
		Scanner refIndScanner = new Scanner(indFile);
		while (refIndScanner.hasNextLine())
		{
			String[] splittedLine = refIndScanner.nextLine().split("\t");
			refIndividualToPopulationMap.put(splittedLine[0], splittedLine[2]);
		}
		refIndScanner.close();
		return refIndividualToPopulationMap;
	}

	static public HashMap<String /*variant*/, HashMap<String /*individual*/, Integer /*eigenstrat genotype*/>> getEigenstratGenotypesFromPlinkPedFile(File pedFile, HashMap<String, String> userIndividualToPopulationMapToFill, LinkedHashSet<Integer> redundantVariantIndexes, String[] variants, HashMap<String, HashMap<String, Integer>> perVariantGtCodes, String progressIndicatorFilePath) throws Exception
	{
		File progressIndicatorFile = new File(progressIndicatorFilePath);
		StringBuffer errors = new StringBuffer();
		
		String sStepLabel = "Reading genotypes provided by user";
		
		LineNumberReader lnr = new LineNumberReader(new FileReader(pedFile));
		lnr.skip(Long.MAX_VALUE);
		int nLineTotalLineCount = 1 + lnr.getLineNumber(), nLineCounter = 0;
		lnr.close();
		
		HashMap<String /*variant*/, HashMap<String /*individual*/, Integer /*eigenstrat genotype*/>> eigenstratGenotypes = new HashMap<String, HashMap<String, Integer>>();
		for (String variant : variants)		// initialize map contents for all variants 
			eigenstratGenotypes.put(variant, new HashMap<String, Integer>());
		HashSet<String> skippedVariants = new HashSet<String>();	// will hold variants for which we have missing data
		
		int nPreviousProgress = 0;
		Scanner pedScanner = new Scanner(pedFile);
		while (pedScanner.hasNextLine())
		{
			String sLine = pedScanner.nextLine().replaceAll("\t", " ");
			if (sLine.trim().length() == 0)
			{
				errors.append("\n- Found empty line in " + pedFile.getName() + " at position " + nLineCounter);
				continue;
			}
			else if (sLine.startsWith("#"))
			{
				LOG.info("Skipping comment at position " + nLineCounter + " in PED file: " + sLine);
				continue;
			}
			
			try
			{
				String sIndividual = readIndividualFromPlinkPedLine(sLine, userIndividualToPopulationMapToFill);
				String[] individualGenotypes = readGenotypesFromPlinkPedLine(sLine, redundantVariantIndexes, variants);
	
				int nCurrentVariantIndex = 0;
				for (String gtCode : individualGenotypes)
				{
					String variant = variants[nCurrentVariantIndex++];
					if (skippedVariants.contains(variant))
						continue;
	
					HashMap<String, Integer> individualToGenotypeMap = eigenstratGenotypes.get(variant);
					HashMap<String /*2-char Plink genotype string*/, Integer /*Eigenstrat-type genotype code*/> eigenstratGtCodes = perVariantGtCodes.get(variant); 
					if (eigenstratGtCodes == null)
						skippedVariants.add(variant);
					else
						individualToGenotypeMap.put(sIndividual, eigenstratGtCodes.get(gtCode));
				}
			}
			catch (WrongNumberArgsException wnae)
			{
				errors.append("\n- " + wnae.getMessage());
				break;
			}

			if (errors.length() > 0)
			{
				pedScanner.close();
				throw new Exception("Uploaded data is invalid: \n" + errors.toString());
			}
			
			if (nLineTotalLineCount > 10)
			{
				int nProgress = (++nLineCounter * 100 / nLineTotalLineCount);
				if (nProgress > nPreviousProgress)
				{
					logProgress(progressIndicatorFile, sStepLabel + "... " + nProgress + "%");
					nPreviousProgress = nProgress;
				}
			}
		}
		pedScanner.close();

		// now remove irrelevant variants from data to return
		for (String aVariantToSkip : skippedVariants)
			eigenstratGenotypes.remove(aVariantToSkip);
		
		if (skippedVariants.size() > 0)
			LOG.info("No data was found for " + (skippedVariants.size()) + " out of " + variants.length + " variants");
		return eigenstratGenotypes;
	}

	static public String readIndividualFromPlinkPedLine(String sPlinkPedLine, Map<String, String> userIndividualToPopulationMapToFill)
	{
		int nFoundSpaceCount = 0, nSpacePos = -1;
		String sIndividual = null, sPopulation = null;
		while (nFoundSpaceCount < 6)
		{
			int nNextSpacePos = sPlinkPedLine.indexOf(" ", nSpacePos + 1);
			if (nFoundSpaceCount == 0)
				sPopulation = sPlinkPedLine.substring(nSpacePos + 1, nNextSpacePos);
			else if (nFoundSpaceCount == 1)
				sIndividual = sPlinkPedLine.substring(nSpacePos + 1, nNextSpacePos);
			nSpacePos = nNextSpacePos;
			nFoundSpaceCount++;
		}
		if (userIndividualToPopulationMapToFill != null)
			userIndividualToPopulationMapToFill.put(sIndividual, sPopulation);
		return sIndividual;
	}
	
	static public String /*3-char string: nucleotide1+space+nucleotide2*/[] readGenotypesFromPlinkPedLine(String sPlinkPedLine, LinkedHashSet<Integer> redundantVariantIndexes, String[] variants) throws WrongNumberArgsException
	{
		String[] result = new String[variants.length];
		
		int nFoundSpaceCount = 0, nSpacePos = -1;
		while (nFoundSpaceCount < 6)
		{
			nSpacePos = sPlinkPedLine.indexOf(" ", nSpacePos + 1);
			nFoundSpaceCount++;
		}
		int nCurrentPos = nSpacePos - 3, nCurrentVariantIndex = -1;
		while (nCurrentPos + 4 < sPlinkPedLine.length())
		{
			nCurrentPos += 4;
			if (!redundantVariantIndexes.contains(nCurrentVariantIndex))	// otherwise it's a duplicate so we ignore it
			{
				nCurrentVariantIndex++;
				if (variants.length + redundantVariantIndexes.size() <= nCurrentVariantIndex)
					throw new WrongNumberArgsException("PED file contains more genotypes than the number of variants defined in MAP file");

				result[nCurrentVariantIndex] = sPlinkPedLine.substring(nCurrentPos, nCurrentPos + 3);
			}
		}
		return result;
	}
		
	static public void logProgress(File f, String s) throws IOException
	{
		LOG.debug(s);
		if (f != null)
			FileUtils.writeStringToFile(f, s);
	}
	
	static public File unzipSingleFileIfNeeded(File zippedFile)	throws Exception
	{
		File result = zippedFile;
		byte[] buffer = new byte[1024];
    	ZipInputStream zis = new ZipInputStream(new FileInputStream(zippedFile));
    	ZipEntry ze = zis.getNextEntry();
    	if (ze != null)
    	{
    	   String fileName = ze.getName();
           result = new File(fileName);
           LOG.debug("file unzip: "+ result.getAbsoluteFile());

           FileOutputStream fos = new FileOutputStream(result);             
           int len;
           while ((len = zis.read(buffer)) > 0)
        	   fos.write(buffer, 0, len);        		
           fos.close();   
    	}
        if (zis.getNextEntry() != null)
        {
            zis.closeEntry();
            zis.close();
        	throw new Exception("Only a single file may be zipped in the genotype file archive!");
        }

        zis.closeEntry();
        zis.close();

        return result;
	}
}
