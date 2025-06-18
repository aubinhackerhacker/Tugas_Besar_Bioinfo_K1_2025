# Tugas Besar Bioinformatika K1 2025
Pipeline dan Script R yang digunakan selama analisis

Pipeline:
![image](https://github.com/user-attachments/assets/a226b63f-f0cd-455a-be8b-967c3eb1b58a)
pipeline ini juga merujuk pada pipeline analisis proteomik (https://github.com/Mr-go-cool/Proteomics-Analysis-Pipeline_1?tab=readme-ov-file), khususnya pada keenam langkah awal

Penjelasan singkat pipeline:
1. Akuisisi data diambil dari database PRIDE dalam bentuk .raw
2. Konversi data dilakukan dari bentuk .raw menjadi bentuk .mzML dengan bantuan ProteoWizard MSConvert (https://proteowizard.sourceforge.io/)
3. Identifikasi protein dilakukan dari file .mzML dengan bantuan DIA-NN (https://github.com/vdemichev/DiaNN) ke dalam bentuk .csv
4. Hasil .csv akan diekstraksi dalam script R yang dapat dilihat pada bagian bawah sehingga didapatkan hasil _gene mapping_ serta _KEGG enrichment analysis_
5. Kedua file akan digunakan untuk analisis lanjutan pada langkah 7 (script R dapat dilihat pada bagian bawah

- Link file yang diunduh untuk dianalisis: https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD053187-1&test=no
- Hasil konversi file .raw ke .mzML: https://drive.google.com/file/d/1oLDz1K8t4GVEupaW7jI3y2i4RyAZ86KE/view?usp=sharing

Repository ini terdiri dari 6 file
1. genemapped.csv (output dari langkah 5)
2. proteins_KEGG_enrichment.csx (output dari langkah 6)
3. Mapped_Gene_Tubes dan KEGG_enrichment_tubes (bentuk .xlsx dari kedua file diatas)
4. Data pre-processing.R (script R yang digunakan untuk pengerjaan langkah 4-6)
5. visualisasi_analisis.R (script R yang digunakan untuk pengerjaan langkah 7)
