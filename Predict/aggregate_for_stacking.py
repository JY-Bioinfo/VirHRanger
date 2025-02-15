import pandas as pd
import argparse
import pickle

'''Global variable'''

'''Order is important!'''
Model_list = {}
Model_list['Genomic'] = ['Genomic_FM', 'Genomic_Traits']
Model_list['Genomic_Protein'] = ['Genomic_FM', 'Protein_FM', 'Genomic_Traits', 'Protein_Traits', 'Homology']
Model_list['Genomic_Protein_PPI'] = ['Genomic_FM', 'Protein_FM', 'Genomic_Traits', 'Protein_Traits', 'Homology','PPI']

LABELS = {}
LABELS['Genomic'] = ['1_Arthropoda', '1_Chordata', '1_Mollusca', '2_Actinopteri',
                 '2_Arachnida', '2_Aves', '2_Insecta', '2_Mammalia', '3_Anseriformes',
                 '3_Artiodactyla', '3_Carnivora', '3_Chiroptera', '3_Diptera',
                 '3_Eulipotyphla', '3_Galliformes', '3_Ixodida', '3_Lagomorpha',
                 '3_Lepidoptera', '3_Passeriformes', '3_Perissodactyla', '3_Primates',
                 '3_Rodentia', '4_Anatidae', '4_Bovidae', '4_Camelidae', '4_Canidae',
                 '4_Cercopithecidae', '4_Cervidae', '4_Cricetidae', '4_Culicidae',
                 '4_Equidae', '4_Felidae', '4_Hominidae', '4_Ixodidae', '4_Leporidae',
                 '4_Muridae', '4_Mustelidae', '4_Phasianidae', '4_Pteropodidae',
                 '4_Rhinolophidae', '4_Sciuridae', '4_Suidae', '4_Vespertilionidae',
                 '5_Aedes', '5_Anas', '5_Bos', '5_Canis', '5_Capra', '5_Chlorocebus',
                 '5_Culex', '5_Equus', '5_Felis', '5_Gallus', '5_Homo', '5_Macaca',
                 '5_Mus', '5_Ovis', '5_Pan', '5_Papio', '5_Pteropus', '5_Rattus',
                 '5_Rhinolophus', '5_Sus', '6_Anas platyrhynchos', '6_Bos taurus',
                 '6_Canis lupus', '6_Capra hircus', '6_Chlorocebus aethiops',
                 '6_Equus caballus', '6_Felis catus', '6_Gallus gallus',
                 '6_Homo sapiens', '6_Macaca fascicularis', '6_Macaca mulatta',
                 '6_Mus musculus', '6_Ovis aries', '6_Pan troglodytes',
                 '6_Rattus norvegicus', '6_Sus scrofa']
LABELS['Genomic_Protein'] = ['1_Arthropoda', '1_Chordata', '1_Mollusca', '2_Actinopteri',
                     '2_Arachnida', '2_Aves', '2_Insecta', '2_Mammalia', '3_Anseriformes',
                     '3_Artiodactyla', '3_Carnivora', '3_Chiroptera', '3_Diptera',
                     '3_Eulipotyphla', '3_Galliformes', '3_Ixodida', '3_Lagomorpha',
                     '3_Lepidoptera', '3_Passeriformes', '3_Perissodactyla', '3_Primates',
                     '3_Rodentia', '4_Anatidae', '4_Bovidae', '4_Camelidae', '4_Canidae',
                     '4_Cercopithecidae', '4_Cervidae', '4_Cricetidae', '4_Culicidae',
                     '4_Equidae', '4_Felidae', '4_Hominidae', '4_Ixodidae', '4_Leporidae',
                     '4_Muridae', '4_Mustelidae', '4_Phasianidae', '4_Pteropodidae',
                     '4_Rhinolophidae', '4_Sciuridae', '4_Suidae', '4_Vespertilionidae',
                     '5_Aedes', '5_Anas', '5_Bos', '5_Canis', '5_Capra', '5_Chlorocebus',
                     '5_Culex', '5_Equus', '5_Felis', '5_Gallus', '5_Homo', '5_Macaca',
                     '5_Mus', '5_Ovis', '5_Pan', '5_Papio', '5_Pteropus', '5_Rattus',
                     '5_Rhinolophus', '5_Sus', '6_Anas platyrhynchos', '6_Bos taurus',
                     '6_Canis lupus', '6_Capra hircus', '6_Chlorocebus aethiops',
                     '6_Equus caballus', '6_Felis catus', '6_Gallus gallus',
                     '6_Homo sapiens', '6_Macaca fascicularis', '6_Macaca mulatta',
                     '6_Mus musculus', '6_Ovis aries', '6_Pan troglodytes',
                     '6_Rattus norvegicus', '6_Sus scrofa']
LABELS['Genomic_Protein_PPI'] = ['6_Homo sapiens']


def parse_args():

    parser = argparse.ArgumentParser(description="prepare data")
    parser.add_argument(
        "--predict_folder",
        type=str,
        required=True,
        help="working folder",
    )
    parser.add_argument(
        "--save_folder",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--ensemble_mode",
        type=str,
        required=True,
    )

    args = parser.parse_args()

    return args


def main():
    args = parse_args()

    PROBS = {}
    six_feat_vids, five_feat_vids, two_feat_vids = [], [], []
    if args.ensemble_mode == 'Genomic_Protein_PPI':
        for MODEL in Model_list[args.ensemble_mode]:
            PROBS[MODEL] = pd.read_csv(
                args.predict_folder + f'{MODEL}/whole_probs.csv')
        '''find intersection for base model'''
        prob_vid_list = []
        for MODEL in Model_list['Genomic_Protein_PPI']:
            prob_vid_list.append(PROBS[MODEL].vid.to_list())
        six_feat_vids = list(set.intersection(*map(set, prob_vid_list)))

        prob_vid_list = []
        for MODEL in Model_list['Genomic_Protein']:
            prob_vid_list.append(PROBS[MODEL].vid.to_list())
        five_feat_vids = list(set.intersection(*map(set, prob_vid_list)))

        prob_vid_list = []
        for MODEL in Model_list['Genomic']:
            prob_vid_list.append(PROBS[MODEL].vid.to_list())
        two_feat_vids = list(set.intersection(*map(set, prob_vid_list)))

        '''subtract'''
        two_feat_vids = [
            item for item in two_feat_vids if item not in five_feat_vids]

        print(
            f'Considering some features are unavailable for some virus genomes, different Ensemble mode are used: \n Genomic_Protein_PPI\t: {len(six_feat_vids)} sequences \n Genomic_Protein\t: {len(five_feat_vids)} sequences \n Genomic\t: {len(two_feat_vids)} sequences \n')

    elif args.ensemble_mode == 'Genomic_Protein':
        for MODEL in Model_list[args.ensemble_mode]:
            PROBS[MODEL] = pd.read_csv(
                args.predict_folder + f'{MODEL}/whole_probs.csv')

        '''find intersection for base model'''
        prob_vid_list = []
        for MODEL in Model_list['Genomic_Protein']:
            prob_vid_list.append(PROBS[MODEL].vid.to_list())
        five_feat_vids = list(set.intersection(*map(set, prob_vid_list)))

        prob_vid_list = []
        for MODEL in Model_list['Genomic']:
            prob_vid_list.append(PROBS[MODEL].vid.to_list())
        two_feat_vids = list(set.intersection(*map(set, prob_vid_list)))

        '''subtract'''
        two_feat_vids = [
            item for item in two_feat_vids if item not in five_feat_vids]

        print(
            f'Considering some features are unavailable for some virus genomes, different Ensemble mode are used: \n Genomic_Protein\t: {len(five_feat_vids)} sequences \n Genomic\t: {len(two_feat_vids)} sequences \n')

    elif args.ensemble_mode == 'Genomic':
        for MODEL in Model_list[args.ensemble_mode]:
            PROBS[MODEL] = pd.read_csv(
                args.predict_folder + f'{MODEL}/whole_probs.csv')

        prob_vid_list = []
        for MODEL in Model_list['Genomic']:
            prob_vid_list.append(PROBS[MODEL].vid.to_list())
        two_feat_vids = list(set.intersection(*map(set, prob_vid_list)))

        print(
            f'Considering some features are unavailable for some virus genomes, different Ensemble mode are used: \n Genomic\t: {len(two_feat_vids)} sequences \n')

    '''Must sort vid!'''
    six_feat_vids = sorted(six_feat_vids)
    five_feat_vids = sorted(five_feat_vids)
    two_feat_vids = sorted(two_feat_vids)

    '''split three modes'''
    six_probs, five_probs, two_probs = {}, {}, {}

    for MODEL in Model_list[args.ensemble_mode]:
        # special process
        six_probs[MODEL] = PROBS[MODEL].loc[PROBS[MODEL].vid.isin(
            six_feat_vids), LABELS['Genomic_Protein_PPI']].reset_index(drop=True)
        if MODEL != 'PPI':
            five_probs[MODEL] = PROBS[MODEL].loc[PROBS[MODEL].vid.isin(
                five_feat_vids), LABELS['Genomic_Protein']].reset_index(drop=True)
            two_probs[MODEL] = PROBS[MODEL].loc[PROBS[MODEL].vid.isin(
                two_feat_vids), LABELS['Genomic']].reset_index(drop=True)

    '''dump files for three modes'''
    if len(six_feat_vids) > 0:
        save_path = args.save_folder + 'Genomic_Protein_PPI/'
        with open(save_path + 'vids.pkl', 'wb') as handle:
            pickle.dump(six_feat_vids, handle)

        with open(save_path + 'test_prob_agg.pkl', 'wb') as handle:
            pickle.dump(six_probs, handle)

        with open(save_path + 'LABELS.pkl', 'wb') as handle:
            pickle.dump(LABELS['Genomic_Protein_PPI'], handle)

        with open(save_path + 'MODEL_list.pkl', 'wb') as handle:
            pickle.dump(Model_list['Genomic_Protein_PPI'], handle)

    if len(five_feat_vids) > 0:
        save_path = args.save_folder + 'Genomic_Protein/'
        with open(save_path + 'vids.pkl', 'wb') as handle:
            pickle.dump(five_feat_vids, handle)

        with open(save_path + 'test_prob_agg.pkl', 'wb') as handle:
            pickle.dump(five_probs, handle)

        with open(save_path + 'LABELS.pkl', 'wb') as handle:
            pickle.dump(LABELS['Genomic_Protein'], handle)

        with open(save_path + 'MODEL_list.pkl', 'wb') as handle:
            pickle.dump(Model_list['Genomic_Protein'], handle)

    if len(two_feat_vids) > 0:
        save_path = args.save_folder + 'Genomic/'
        with open(save_path + 'vids.pkl', 'wb') as handle:
            pickle.dump(two_feat_vids, handle)

        with open(save_path + 'test_prob_agg.pkl', 'wb') as handle:
            pickle.dump(two_probs, handle)

        with open(save_path + 'LABELS.pkl', 'wb') as handle:
            pickle.dump(LABELS['Genomic'], handle)

        with open(save_path + 'MODEL_list.pkl', 'wb') as handle:
            pickle.dump(Model_list['Genomic'], handle)


if __name__ == "__main__":
    main()
