# **Notes**

## To-do:
1. Make new SRA submission for BioProject PRJNA693865 adding mock + NTC fastq files and associated metadata.
2. Prepare public version of FS9 repo saving only the relevant code (all old code - save in separate repo)

## BioProject PRJNA693865 SRA submission modification

### 25Jan2021: Submitted SRA submission SUB8884056

### 1Mar2021
* Realized I included DNEG3 samples in SRA submission SUB8884056.
* Emailed NLM Support how to make changes, including:
  * BioSample attributes file (FS9_Model.organism.animal.1.0.xlsx)
  * Metadata file (FS9_SRA_metadata.xlsx)
  * Removal (DNEG3) and addition (mock, NTC) of several fastq files.

```
Dear Kathy ,
1.  BioSample attributes file - please find the file in attachment, edit it and send back to me
2. you can make spreadsheet-like updates of SRA metadata  in DataView  interface  https://dataview.ncbi.nlm.nih.gov/object/PRJNA693865
3. removal and addition of several fastq files -  we cannot replace loaded files, please provide run accessions that you want to be withdrawn; to add new sequence data to the bioproject (and/or biosamples)  please make new SRA submission and use this bioproject's and corresponding biosamples' (if applicable) accessions for SRA metadata. Please note that you either can use existing biosamples in the SRA submission or create all new ones - you cannot mix and match.

For more information please refer to the Update guide: https://www.ncbi.nlm.nih.gov/sra/docs/submitupdate/

Best regards,
Svetlana

If you have any questions or concerns regarding your SRA submission please don’t hesitate to contact sra@ncbi.nlm.nih.gov (applies to new questions).  We normally respond within 2 business days.

Svetlana Iazvovskaia
The NCBI SRA database submission staff
```

### 9Mar2021
Follow-up response from Svetlana when I was still unclear what I could and couldn't do with the SRA submission (needed clarity):
```
Hi Kathy,

You cannot submit same SRA data again - i.e. you cannot withdraw SRA from one submission and then, submit them again in another.
If you want to keep SRA data but do not like how samples look, please create new samples and relink your SRA data to new samples and let me know when you done I will withdraw old samples (that will be unlinked from SRA).

If you want to withdraw SRA data (permanently - not to ever submit them again) - please let me know SRRs, if you'd like to withdraw the linked samples as well - please let me know.

If you want to add new SRA and new samples - please create new SRA submission.

Please let me know what you'd like to do.

Not knowing nature of your samples' modifications I would avoid editing old samples and prefer to withdraw them and for you to create new ones.  Samples are objects. Minor change  to an object is ok but converting one object to completely different object is not.

Hope it helps.

Please let me know what you want to withdraw - runs (SRR) and samples (SAMNs) - when you are ready.

Best regards,

Svetlana
```

### 17Mar2021
Received this message. Afterwards, submitted modified BioSample attributes file `Copy of biosamples.195.2021-03-01T14-41-30.xlsx` with corrections to days. Need to upload new SRA submission of mock and NTC samples.
```
Hi Kathy,

I had withdrawn samples and SRA in the list.
Please send me list of samples with the modified attribute.
Please find spreadsheet for remaining samples in the attachment..

Best regards,

Svetlana
```
