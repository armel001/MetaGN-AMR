rgi clean --local
rgi load \
  --card_json /data/armel/git/MetagenAMR-main/resources/databases/card_db/card.json \
  --debug --local \
  --card_annotation /data/armel/git/MetagenAMR-main/resources/databases/card_db/card_database_v4.0.1.fasta \
  --card_annotation_all_models /data/armel/git/MetagenAMR-main/resources/databases/card_db/card_database_v4.0.1_all.fasta \
  --wildcard_annotation /data/armel/git/MetagenAMR-main/resources/databases/card_db/wildcard_database_v4.0.1.fasta \
  --wildcard_annotation_all_models /data/armel/git/MetagenAMR-main/resources/databases/card_db/wildcard_database_v4.0.1_all.fasta \
  --wildcard_index /data/armel/git/MetagenAMR-main/resources/databases/card_db/wildcard/index-for-model-sequences.txt \
  --wildcard_version 4.0.1 \
  --amr_kmers /data/armel/git/MetagenAMR-main/resources/databases/card_db/wildcard/all_amr_61mers.txt \
  --kmer_database /data/armel/git/MetagenAMR-main/resources/databases/card_db/wildcard/61_kmer_db.json \
  --kmer_size 61
