import torch
import torch.utils.checkpoint
from torch import nn, Tensor
import torch.nn.functional as F
from torch.nn import BCEWithLogitsLoss
import math


class MHC_Net_Head(nn.Module):
    def __init__(self, config):

        super().__init__()
        self.hierarchical_depth = config.hierarchical_depth
        self.hierarchical_class = config.hierarchy_classes
        self.g_size = config.g_size
        self.G2L = config.G2L
        self.hidden_size = config.hidden_size
        self.local_layers = torch.nn.ModuleList()
        self.global_layers = torch.nn.ModuleList()
        self.num_labels = config.num_labels
        self.cls_dropout_prob = config.cls_dropout_prob

        self.global1 = torch.nn.Linear(self.hidden_size, self.g_size)
        self.global2 = torch.nn.Linear(
            self.g_size + self.hidden_size, self.g_size)
        self.global3 = torch.nn.Linear(
            self.g_size + self.hidden_size, self.g_size)
        self.global4 = torch.nn.Linear(
            self.g_size + self.hidden_size, self.g_size)
        self.global5 = torch.nn.Linear(
            self.g_size + self.hidden_size, self.g_size)
        self.global6 = torch.nn.Linear(
            self.g_size + self.hidden_size, self.g_size)

        self.out_proj_G = torch.nn.Linear(self.g_size, self.num_labels)

        self.local1 = torch.nn.Linear(self.g_size, self.G2L)
        self.local2 = torch.nn.Linear(self.g_size, self.G2L)
        self.local3 = torch.nn.Linear(self.g_size, self.G2L)
        self.local4 = torch.nn.Linear(self.g_size, self.G2L)
        self.local5 = torch.nn.Linear(self.g_size, self.G2L)
        self.local6 = torch.nn.Linear(self.g_size, self.G2L)

        self.out_proj_L1 = torch.nn.Linear(
            self.G2L, self.hierarchical_class[0])
        self.out_proj_L2 = torch.nn.Linear(
            self.G2L, self.hierarchical_class[1])
        self.out_proj_L3 = torch.nn.Linear(
            self.G2L, self.hierarchical_class[2])
        self.out_proj_L4 = torch.nn.Linear(
            self.G2L, self.hierarchical_class[3])
        self.out_proj_L5 = torch.nn.Linear(
            self.G2L, self.hierarchical_class[4])
        self.out_proj_L6 = torch.nn.Linear(
            self.G2L, self.hierarchical_class[5])

        self.dropout = torch.nn.Dropout(self.cls_dropout_prob)

    def forward(self, features, **kwargs):
        x = features

        # start gradient flow

        G1 = self.dropout(F.relu(self.global1(x)))
        L1 = self.dropout(F.relu(self.local1(G1)))
        L1_Out = self.out_proj_L1(L1)

        G2 = self.dropout(F.relu(self.global2(torch.cat([G1, x], dim=1))))
        L2 = self.dropout(F.relu(self.local2(G2)))
        L2_Out = self.out_proj_L2(L2)

        G3 = self.dropout(F.relu(self.global3(torch.cat([G2, x], dim=1))))
        L3 = self.dropout(F.relu(self.local3(G3)))
        L3_Out = self.out_proj_L3(L3)

        G4 = self.dropout(F.relu(self.global4(torch.cat([G3, x], dim=1))))
        L4 = self.dropout(F.relu(self.local4(G4)))
        L4_Out = self.out_proj_L4(L4)

        G5 = self.dropout(F.relu(self.global5(torch.cat([G4, x], dim=1))))
        L5 = self.dropout(F.relu(self.local5(G5)))
        L5_Out = self.out_proj_L5(L5)

        G6 = self.dropout(F.relu(self.global6(torch.cat([G5, x], dim=1))))
        L6 = self.dropout(F.relu(self.local6(G6)))
        L6_Out = self.out_proj_L6(L6)

        G_Out = self.out_proj_G(G6)
        L_Out = torch.cat(
            [L1_Out, L2_Out, L3_Out, L4_Out, L5_Out, L6_Out], dim=1)

        return G_Out, L_Out, 0.5*G_Out + 0.5*L_Out


def violation_pen(hierarchy_structure, weighted_logits):
    probs = F.softmax(weighted_logits, dim=1)
    pens = []
    for child in hierarchy_structure.keys():
        parent = hierarchy_structure[child]
        p_child = probs[:, child]
        p_parent = probs[:, parent]
        temp = F.relu(p_child-p_parent)
        # mean for the whole batch
        pens.append(torch.mean(torch.pow(input=temp, exponent=2)))
    # sum for all the labels
    return torch.sum(torch.stack(pens))


class MHC_Net(nn.Module):

    def __init__(self, config):

        super().__init__()
        self.classifier = MHC_Net_Head(config)
        # input projection
        self.in_proj = nn.Linear(config.feature_dim, config.hidden_size)
        self.num_labels = config.num_labels
        # add for penalty loss
        self.loss_type = config.loss_type
        self.criterion = BCEWithLogitsLoss
        self.hierarchy_structure = config.hierarchy_structure
        self.coe_l = config.coe_l

    def forward(
        self,
        features=None,
        labels=None,
        weights=None,
    ):
        proj_features = self.in_proj(features)

        global_logits, local_logits, weighted_logits = self.classifier(
            proj_features)

        if labels is None:
            loss = 0
        else:

            if self.loss_type == 'Focal':
                loss_fct = Binary_FocalLoss(weights=weights)
                g_loss = loss_fct(global_logits, labels.float())
                l_loss = loss_fct(local_logits, labels.float())
                loss = g_loss + l_loss

            elif self.loss_type == 'penalty':
                loss_fct = self.criterion(weight=weights)
                g_loss = loss_fct(global_logits, labels.float())
                l_loss = loss_fct(local_logits, labels.float())
                p_loss = violation_pen(
                    self.hierarchy_structure, weighted_logits)
                loss = g_loss + l_loss + self.coe_l * p_loss

            else:
                loss_fct = self.criterion(weight=weights)
                g_loss = loss_fct(global_logits.float(), labels.float())
                l_loss = loss_fct(local_logits.float(), labels.float())
                # avoid nan for logits
                if g_loss != g_loss:
                    print(global_logits.float())

                loss = g_loss + l_loss

        return_dict = {}
        return_dict['logits'] = weighted_logits
        return_dict['loss'] = loss

        return return_dict


class FlatHeadSix(nn.Module):

    """Head for sentence-level classification tasks."""

    def __init__(self, config):
        super().__init__()
        self.num_labels = config.num_labels
        self.g_size = config.g_size
        self.hidden_size = config.hidden_size
        self.cls_dropout_prob = config.cls_dropout_prob
        self.dropout = torch.nn.Dropout(self.cls_dropout_prob)

        self.D1 = nn.Linear(self.hidden_size, self.g_size)
        self.D2 = nn.Linear(self.g_size, self.g_size)
        self.D3 = nn.Linear(self.g_size, self.g_size)
        self.D4 = nn.Linear(self.g_size, self.g_size)
        self.D5 = nn.Linear(self.g_size, self.g_size)
        self.D6 = nn.Linear(self.g_size, self.g_size)
        self.out_proj = nn.Linear(self.g_size, self.num_labels)

    def forward(self, x, **kwargs):

        x = self.dropout(F.relu(self.D1(x)))
        x = self.dropout(F.relu(self.D2(x)))
        x = self.dropout(F.relu(self.D3(x)))
        x = self.dropout(F.relu(self.D4(x)))
        x = self.dropout(F.relu(self.D5(x)))
        x = self.dropout(F.relu(self.D6(x)))
        x = self.out_proj(x)

        return x


class Flat_Net(nn.Module):

    def __init__(self, config):

        super().__init__()
        self.classifier = FlatHeadSix(config)
        # input projection
        self.in_proj = nn.Linear(config.feature_dim, config.hidden_size)
        self.num_labels = config.num_labels
        self.loss_type = config.loss_type
        self.criterion = BCEWithLogitsLoss

    def forward(
        self,
        features=None,
        labels=None,
        weights=None,
    ):
        proj_features = self.in_proj(features)

        logits = self.classifier(proj_features)

        loss = None
        if labels is None:
            loss = 0

        else:
            if self.loss_type == 'Focal':
                loss_fct = Binary_FocalLoss(weights=weights)
                loss = loss_fct(logits.float(), labels.float())

            else:
                loss_fct = self.criterion(weight=weights)
                loss = loss_fct(logits.float(), labels.float())

        return_dict = {}
        return_dict['logits'] = logits
        return_dict['loss'] = loss

        return return_dict


class PositionalEncoding(nn.Module):

    def __init__(self, in_dim, max_len, config):

        super().__init__()
        self.dim = in_dim
        self.dropout = nn.Dropout(config.dropout)
        self.max_len = max_len

        position = torch.arange(self.max_len).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, self.dim, 2)
                             * (-math.log(10000.0) / self.dim))
        pe = torch.zeros(self.max_len, 1, self.dim)
        pe[:, 0, 0::2] = torch.sin(position * div_term)
        pe[:, 0, 1::2] = torch.cos(position * div_term)
        self.register_buffer('pe', pe)

    def forward(self, x: Tensor) -> Tensor:
        """
        Args:
            x: Tensor, shape [seq_len, batch_size, embedding_dim]
        """
        x = x + self.pe[:x.size(0)]
        return self.dropout(x)


class TransformerBlock(nn.Module):
    def __init__(self, in_dim, max_len, config):
        super().__init__()
        self.model_type = 'Transformer'
        self.dim = in_dim
        self.max_len = max_len
        self.heads = config.heads
        self.dropout = config.dropout
        self.num_encoder_layers = config.num_encoder_layers

        self.pos_embedding = PositionalEncoding(self.dim, self.max_len, config)
        encoder_layer = nn.TransformerEncoderLayer(
            self.dim, self.heads, dropout=self.dropout)
        encoder_norm = nn.LayerNorm(self.dim)
        self.encoder = nn.TransformerEncoder(
            encoder_layer, self.num_encoder_layers, encoder_norm)

    def generate_square_subsequent_mask(self, sz):
        mask = (torch.triu(torch.ones(sz, sz)) == 1).transpose(0, 1)
        mask = mask.float().masked_fill(mask == 0, float(
            '-inf')).masked_fill(mask == 1, float(0.0))
        return mask  # [sz,sz]

    def forward(self, src,
                src_key_padding_mask=None,
                concat_type='sum'):

        src_embed = self.pos_embedding(src)  # [src_len,batch_size,embed_dim]
        # [src_len,batch_size,embed_dim]
        memory = self.encoder(
            src=src_embed, src_key_padding_mask=src_key_padding_mask)
        # print(memory.shape)

        # pooling on the sequence L dimension
        if concat_type == 'sum':
            memory = torch.sum(memory, dim=0)
        elif concat_type == 'avg':
            memory = torch.sum(memory, dim=0) / memory.size(0)
        else:
            memory = memory[-1, ::]

        return memory


class MHC_hybrid_Net(nn.Module):
    def __init__(self, config):
        super().__init__()

        self.model_type = 'Hybrid'
        self.nc_dim = config.nc_dim
        self.pro_dim = config.pro_dim
        self.dropout = config.dropout
        self.num_labels = config.num_labels

        self.nc_encoder = TransformerBlock(
            self.nc_dim, config.nc_max_len, config)
        self.pro_encoder = TransformerBlock(
            self.pro_dim, config.pro_max_len, config)
        self.classifier = MHC_Net_Head(config)

        # add for penalty loss
        self.loss_type = config.loss_type
        self.criterion = BCEWithLogitsLoss
        self.hierarchy_structure = config.hierarchy_structure
        self.coe_l = config.coe_l

    def _init_weights(self, module):
        """ Initialize the weights """
        if isinstance(module, nn.Linear):
            module.weight.data.normal_(mean=0.0, std=0.02)
            if module.bias is not None:
                module.bias.data.zero_()
        elif isinstance(module, nn.LayerNorm):
            module.bias.data.zero_()
            module.weight.data.fill_(1.0)

    def forward(self, nc_feature=None,
                nc_padding_mask=None,
                pro_feature=None,
                pro_padding_mask=None,
                labels=None,
                weights=None,
                concat_type='sum'):

        nc_feature = nc_feature.transpose(0, 1)
        pro_feature = pro_feature.transpose(0, 1)
        if nc_padding_mask != None:
            nc_padding_mask = nc_padding_mask.bool()
        if pro_padding_mask != None:
            pro_padding_mask = pro_padding_mask.bool()

        # feature [src_len,batch_size,embed_dim]
        # padding_mask [batch_size, src_len]

        # genome

        nc_memory = self.nc_encoder(nc_feature, nc_padding_mask)
        # need transpose
        # [batch_size,embed_dim]

        # protein
        pro_memory = self.pro_encoder(pro_feature, pro_padding_mask)
        # [batch_size,embed_dim]

        memory = torch.cat((nc_memory, pro_memory), 1)

        global_logits, local_logits, weighted_logits = self.classifier(memory)

        if self.loss_type == 'Focal':
            loss_fct = Binary_FocalLoss(weights=weights)
            g_loss = loss_fct(global_logits, labels.float())
            l_loss = loss_fct(local_logits, labels.float())
            loss = g_loss + l_loss

        elif self.loss_type == 'penalty':
            loss_fct = self.criterion(weight=weights)
            g_loss = loss_fct(global_logits, labels.float())
            l_loss = loss_fct(local_logits, labels.float())
            p_loss = violation_pen(self.hierarchy_structure, weighted_logits)
            loss = g_loss + l_loss + self.coe_l * p_loss

        else:
            loss_fct = self.criterion(weight=weights)
            g_loss = loss_fct(global_logits.float(), labels.float())
            l_loss = loss_fct(local_logits.float(), labels.float())
            # avoid nan for logits
            if g_loss != g_loss:
                print(global_logits.float())

            loss = g_loss + l_loss

        return_dict = {}
        return_dict['logits'] = weighted_logits
        return_dict['loss'] = loss

        return return_dict

        
class Binary_FocalLoss(nn.Module):
    def __init__(self, weights=None, gamma=2, alpha=0.25):
        super(Binary_FocalLoss, self).__init__()
        self.weights = weights
        self.gamma = gamma  # Focusing parameter (γ)
        self.alpha = alpha  # Balancing factor for class imbalance

    def forward(self, inputs, targets):
        # Apply BCEWithLogitsLoss to compute the base loss
        bce_loss = F.binary_cross_entropy_with_logits(inputs, targets, weight=self.weights, reduction='none')
        
        # Apply the focal loss modulating factor
        sigmoid_p = torch.sigmoid(inputs)  # Predicted probabilities
        p_t = targets * sigmoid_p + (1 - targets) * (1 - sigmoid_p)  # Probability of true class
        focal_loss = self.alpha * (1 - p_t) ** self.gamma * bce_loss
        
        return focal_loss.mean()  # Return the mean loss over the batch