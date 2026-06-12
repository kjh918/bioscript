import os 
import sys 
from glob import glob
import pandas as pd 
from pathlib import Path
import argparse
import gzip

def read_recept_data(recept_file_path):
    recept_df = pd.read_csv(recept_file_path, sep='\t')
    if 'Order_no' in recept_df.columns:
        recept_df['Order_no'] = recept_df['Order_no'].astype(str).str.zfill(2)
    return recept_df

def parse_flowcell_id_from_fastq(fastq_path):
    try:
        with gzip.open(fastq_path, 'rt') as f:
            header = f.readline().strip()
            parts = header.split(':')
            if len(parts) >= 3:
                return parts[2]
    except Exception as e:
        print(f"  [Warning] Failed to parse flowcell ID from {fastq_path}: {e}")
    return ""

def link_data(raw_data_path, recept_df, batch_id, workdir_path):
    if 'Data' in recept_df.columns:
        received_data_df = recept_df[recept_df['Data'] == 'O']
    else:
        received_data_df = recept_df 
    
    target_dir = os.path.join(workdir_path, 'FASTQ')
    os.makedirs(target_dir, exist_ok=True)

    for _, row in received_data_df.iterrows():
        data_id = str(row.get('data_id', ''))
        order_num = str(row.get('Order_no', ''))
        
        if not data_id or not order_num:
            continue

        link_prefix = os.path.join(target_dir, f'{batch_id}_{order_num}')
        file_list = glob(os.path.join(raw_data_path, f'{data_id}*fastq.gz'))

        for path in file_list:
            if '_1' in path or '_R1' in path:
                suffix = '_R1.fastq.gz'
            else:
                suffix = '_R2.fastq.gz'
                
            link_path = f"{link_prefix}{suffix}"
            
            if not os.path.exists(link_path):
                cmd = f'ln -s {path} {link_path}'
                os.system(cmd)

def preprocess_metadata(recept_df, batch_id, raw_data_path):
    analysis_df_dict = {
        'seq_id':[], 'panel':[], 'depth':[], 'seq_folder':[], 'flowcell_id':[], 'sequencer_id':[], 'sample_no':[],
        'project_facility_id':[], 'ngs_facility_id':[], 'client_facility_id':[],
        'sample_name':[], 'sample_label':[], 'sample_group':[], 'patient_id':[],
        'sample_type':[], 'sample_origin':[], 'sample_tissue':[], 'sample_info':[],
        'matched_normal':[], 'date_sample_in':[], 'date_sample_qc':[], 'date_ngs_data':[], 'year':[]
    }

    for _, row in recept_df.iterrows():
        order_no = str(row.get('Order_no', ''))
        data_id = str(row.get('data_id', ''))
        
        # --- Flowcell ID 추출 ---
        flowcell_id = ""
        if raw_data_path and data_id:
            file_list = glob(os.path.join(raw_data_path, f'{data_id}*fastq.gz'))
            if file_list:
                flowcell_id = parse_flowcell_id_from_fastq(file_list[0])
        if not flowcell_id:
            flowcell_id = ""

        # =====================================================================
        # 1. 이름 파싱 (recept.tsv에서 patient_id 직접 분리)
        # =====================================================================
        # sample_label: 'Sample Name' 컬럼 값 전체 (예: 09-0002 (RCC))
        sample_label = str(row.get('Sample Name', '')).strip()
        
        # [수정] patient_id: TSV의 'Patient ID' 컬럼에서 직접 가져옴
        patient_id = str(row.get('Patient ID', '')).strip()
        
        # (안전 장치) 만약 TSV에 Patient ID를 비워뒀다면 라벨에서 자동 추출
        if not patient_id:
            if '(' in sample_label:
                patient_id = sample_label.split('(')[0].strip()
            elif len(sample_label) >= 9 and sample_label.startswith('102'):
                patient_id = sample_label[:9]
            else:
                patient_id = sample_label
            
        # sample_group은 patient_id를 따라가도록 설정
        sample_group = patient_id

        # 2. Tissue & Cancer Info, Specimen 파싱
        cancer_type = str(row.get('Cancer Type', '')).strip()
        sample_info = f"h{cancer_type}" if cancer_type else ""
        sample_tissue = 'kidney' if cancer_type == 'RCC' else ''
        sample_origin = str(row.get('Specimen', '')).strip()

        # sample_name 최종 조립: hRCC_09-0002_ORG 형태 (patient_id 기준)
        sample_name = f"{sample_info}_{patient_id}_{sample_origin}"
        # =====================================================================

        # 3. 날짜 변환 (YYMMDD -> YYYYMMDD)
        reception_date_raw = str(row.get('Reception Date', '')).split('.')[0].strip()
        if len(reception_date_raw) == 6:
            date_sample_in = '20' + reception_date_raw
        else:
            date_sample_in = reception_date_raw
            
        year = date_sample_in[:4] if len(date_sample_in) >= 4 else '2026'

        # Dictionary Append
        analysis_df_dict['seq_id'].append(f"{batch_id}_{order_no}")
        analysis_df_dict['panel'].append('WES')
        analysis_df_dict['depth'].append(100 if str(row.get('Depth/Ouput', '')).strip() == '6' else 200)
        analysis_df_dict['seq_folder'].append(batch_id)
        
        analysis_df_dict['flowcell_id'].append(flowcell_id)
        analysis_df_dict['sequencer_id'].append(3)
        analysis_df_dict['sample_no'].append(int(order_no))
        analysis_df_dict['project_facility_id'].append(1)
        analysis_df_dict['ngs_facility_id'].append(14)
        analysis_df_dict['client_facility_id'].append(row.get('Order Facility', ''))
        
        # 파싱된 식별자 매핑
        analysis_df_dict['sample_name'].append(sample_name)
        analysis_df_dict['sample_label'].append(sample_label) 
        analysis_df_dict['sample_group'].append(sample_group) 
        analysis_df_dict['patient_id'].append(patient_id)
        
        analysis_df_dict['sample_type'].append('gDNA') # Cell -> gDNA 고정
        analysis_df_dict['sample_origin'].append(sample_origin)
        analysis_df_dict['sample_tissue'].append(sample_tissue)
        analysis_df_dict['sample_info'].append(sample_info)
        analysis_df_dict['matched_normal'].append('0')
        
        analysis_df_dict['date_sample_in'].append(date_sample_in)
        analysis_df_dict['date_sample_qc'].append('20260409') 
        analysis_df_dict['date_ngs_data'].append('20260512')  
        analysis_df_dict['year'].append(year)

    meta_df = pd.DataFrame(analysis_df_dict)
    return meta_df

def main():
    parser = argparse.ArgumentParser(description="Run organoid pipeline metadata preprocessing")
    
    parser.add_argument('--input', '-i', required=True, help='input file path (recept.tsv)')
    parser.add_argument('--batch_id', '-b', required=True, help='Batch ID for the run (e.g., WES_26_01)')
    parser.add_argument('--rawdata_path', '-r', help='raw data directory path')
    parser.add_argument('--script_path', '-s', default='/storage/home/jhkim/Projects/NGS/WES/v1', help='script directory path')
    parser.add_argument('--workdir_path', '-w', default='/data/wes', help='work directory path')
    parser.add_argument('--report_path', '-R', default='/storage/home/kangsm/shinyWeb/REPORT_WES', help='report file path')

    args = parser.parse_args()
    workdir_path = Path(args.workdir_path)
    
    print(f"[INFO] Reading input file: {args.input}")
    recept_df = read_recept_data(args.input)
    
    print(f"[INFO] Generating metadata for Batch: {args.batch_id}")
    meta_df = preprocess_metadata(recept_df, args.batch_id, args.rawdata_path)
    
    os.makedirs(workdir_path, exist_ok=True)
    output_meta_file = workdir_path / f"{args.batch_id}_metadata.tsv"
    meta_df.to_csv(output_meta_file, sep='\t', index=False, encoding='utf-8')
    print(f"[INFO] Metadata successfully saved to: {output_meta_file}")
    
    if args.rawdata_path:
        print(f"[INFO] Linking raw data from: {args.rawdata_path}")
        link_data(args.rawdata_path, recept_df, args.batch_id, str(workdir_path))

if __name__ == "__main__":
    main()