from deconstruct_lc import tools_lc
from deconstruct_lc.scores.norm_score import NormScore

def main():
    k = 6
    lce = 1.6
    lca = 'SGEQAPDTNKR'
    sup35 = "MSDSNQGNNQQNYQQYSQNGNQQQGNNRYQGYQAYNAQAQPAGGYYQNYQGYSGYQQGGY" \
            "QQYNPDAGYQQQYNPQGGYQQYNPQGGYQQQFNPQGGRGNYKNFNYNNNLQGYQAGFQPQ" \
            "SQGMSLNDFQKQQKQAAPKPKKTLKLVSSSGIKLANATKKVGTKPAESDKKEEEKSAETK" \
            "EPTKEPTKVEEPVKKEEKPVQTEEKTEEKSELPKVEDLKISESTHNTNNANVTSADALIK" \
            "EQEEEVDDEVVNDMFGGKDHVSLIFMGHVDAGKSTMGGNLLYLTGSVDKRTIEKYEREAK" \
            "DAGRQGWYLSWVMDTNKEERNDGKTIEVGKAYFETEKRRYTILDAPGHKMYVSEMIGGAS" \
            "QADVGVLVISARKGEYETGFERGGQTREHALLAKTQGVNKMVVVVNKMDDPTVNWSKERY" \
            "DQCVSNVSNFLRAIGYNIKTDVVFMPVSGYSGANLKDHVDPKECPWYTGPTLLEYLDTMN" \
            "HVDRHINAPFMLPIAAKMKDLGTIVEGKIESGHIKKGQSTLLMPNKTAVEIQNIYNETEN" \
            "EVDMAMCGEQVKLRIKGVEEEDISPGFVLTSPKNPIKSVTKFVAQIAIVELKSIIAAGFS" \
            "CVMHVHTAIEEVHIVKLLHKLEKGTNRKSKKPPAFAKKGMKVIAVLETEAPVCVETYQDY" \
            "PQLGRFTLRDQGTTIAIGKIVKIAE"
    #print(sup35[:260])
    ns = tools_lc.display_lc(sup35[:260], k, lca, lce)
    #print(ns)
    motifs = tools_lc.count_lc_motifs(sup35[:260], k, lca, lce)
    #print(motifs)
    #norm = NormScore([sup35[:260]])
    #print(norm.lc_norm_score())

    pab1 = 'MADITDKTAEQLENLNIQDDQKQAATGSESQSVENSSASLYVGDLEPSVSEAHLYDIFSP' \
           'IGSVSSIRVCRDAITKTSLGYAYVNFNDHEAGRKAIEQLNYTPIKGRLCRIMWSQRDPSL' \
           'RKKGSGNIFIKNLHPDIDNKALYDTFSVFGDILSSKIATDENGKSKGFGFVHFEEEGAAK' \
           'EAIDALNGMLLNGQEIYVAPHLSRKERDSQLEETKAHYTNLYVKNINSETTDEQFQELFA' \
           'KFGPIVSASLEKDADGKLKGFGFVNYEKHEDAVKAVEALNDSELNGEKLYVGRAQKKNER' \
           'MHVLKKQYEAYRLEKMAKYQGVNLFVKNLDDSVDDEKLEEEFAPYGTITSAKVMRTENGK' \
           'SKGFGFVCFSTPEEATKAITEKNQQIVAGKPLYVAIAQRKDVRRSQLAQQIQARNQMRYQ' \
           'QATAAAAAAAAGMPGQFMPPMFYGVMPPRGVPFNGPNPQQMNPMGGMPKNGMPPQFRNGP' \
           'VYGVPPQGGFPRNANDNNQFYQQKQRQALGEQLYKKVSAKTSNEEAAGKITGMILDLPPQ' \
           'EVFPLLESDELFEQHYKEASAAYESFKKEQEQQTEQA'
    wop = pab1[0:420] + pab1[505:]
    print(wop)
    #print(wop)
    motifs = tools_lc.count_lc_motifs(wop, k, lca, lce)
    #print(motifs)
    norm = NormScore()
    print(norm.lc_norm_score([wop]))
    ns = tools_lc.display_lc(wop, k, lca, lce)
    print(ns)
    nck = 'MAEEVVVVAKFDYVAQQEQELDIKKNERLWLLDDSKSWWRVRNSMNKTGFVPSNYVERKN' \
          'SARKASIVKNLKDTLGIGKVKRKPSVPDSASPADDSFVDPGERLYDLNMPAYVKFNYMAE' \
          'REDELSLIKGTKVIVMEKCSDGWWRGSYNGQVGWFPSNYVTEEGDSPLGDHVGSLSEKLA' \
          'AVVNNLNTGQVLHVVQALYPFSSSNDEELNFEKGDVMDVIEKPENDPEWWKCRKINGMVG' \
          'LVPKNYVTVMQNNPLTSGLEPSPPQCDYIRPSLTGKFAGNPWYYGKVTRHQAEMALNERG' \
          'HEGDFLIRDSESSPNDFSVSLKAQGKNKHFKVQLKETVYCIGQRKFSTMEELVEHYKKAP' \
          'IFTSEQGEKLYLVKHLS'
    #norm = NormScore([nck])
    #print(norm.lc_norm_score())
    nwasp = 'MSSVQQQPPPPRRVTNVGSLLLTPQENESLFTFLGKKCVTMSSAVVQLYAADRNCMWSKK' \
            'CSGVACLVKDNPQRSYFLRIFDIKDGKLLWEQELYNNFVYNSPRGYFHTFAGDTCQVALN' \
            'FANEEEAKKFRKAVTDLLGRRQRKSEKRRDPPNGPNLPMATVDIKNPEITTNRFYGPQVN' \
            'NISHTKEKKKGKAKKKRLTKADIGTPSNFQHIGHVGWDPNTGFDLNNLDPELKNLFDMCG' \
            'ISEAQLKDRETSKVIYDFIEKTGGVEAVKNELRRQAPPPPPPSRGGPPPPPPPPHNSGPP' \
            'PPPARGRGAPPPPPSRAPTAAPPPPPPSRPSVAVPPPPPNRMYPPPPPALPSSAPSGPPP' \
            'PPPSVLGVGPVAPPPPPPPPPPPGPPPPPGLPSDGDHQVPTTAGNKAALLDQIREGAQLK' \
            'KVEQNSRPVSCSGRDALLDQIRQGIQLKSVADGQESTPPTPAPTSGIVGALMEVMQKRSK' \
            'AIHSSDEDEDEDDEEDFEDDDEWED'
    #norm = NormScore([nwasp])
    #print(norm.lc_norm_score())
    nephrin = 'MALGTTLRASLLLLGLLTEGLAQLAIPASVPRGFWALPENLTVVEGASVELRCGVSTPGS' \
              'AVQWAKDGLLLGPDPRIPGFPRYRLEGDPARGEFHLHIEACDLSDDAEYECQVGRSEMGP' \
              'ELVSPRVILSILVPPKLLLLTPEAGTMVTWVAGQEYVVNCVSGDAKPAPDITILLSGQTI' \
              'SDISANVNEGSQQKLFTVEATARVTPRSSDNRQLLVCEASSPALEAPIKASFTVNVLFPP' \
              'GPPVIEWPGLDEGHVRAGQSLELPCVARGGNPLATLQWLKNGQPVSTAWGTEHTQAVARS' \
              'VLVMTVRPEDHGAQLSCEAHNSVSAGTQEHGITLQVTFPPSAIIILGSASQTENKNVTLS' \
              'CVSKSSRPRVLLRWWLGWRQLLPMEETVMDGLHGGHISMSNLTFLARREDNGLTLTCEAF' \
              'SEAFTKETFKKSLILNVKYPAQKLWIEGPPEGQKLRAGTRVRLVCLAIGGNPEPSLMWYK' \
              'DSRTVTESRLPQESRRVHLGSVEKSGSTFSRELVLVTGPSDNQAKFTCKAGQLSASTQLA' \
              'VQFPPTNVTILANASALRPGDALNLTCVSVSSNPPVNLSWDKEGERLEGVAAPPRRAPFK' \
              'GSAAARSVLLQVSSRDHGQRVTCRAHSAELRETVSSFYRLNVLYRPEFLGEQVLVVTAVE' \
              'QGEALLPVSVSANPAPEAFNWTFRGYRLSPAGGPRHRILSSGALHLWNVTRADDGLYQLH' \
              'CQNSEGTAEARLRLDVHYAPTIRALQDPTEVNVGGSVDIVCTVDANPILPGMFNWERLGE' \
              'DEEDQSLDDMEKISRGPTGRLRIHHAKLAQAGAYQCIVDNGVAPPARRLLRLVVRFAPQV' \
              'EHPTPLTKVAAAGDSTSSATLHCRARGVPNIVFTWTKNGVPLDLQDPRYTEHTYHQGGVH' \
              'SSLLTIANVSAAQDYALFTCTATNALGSDQTNIQLVSISRPDPPSGLKVVSLTPHSVGLE' \
              'WKPGFDGGLPQRFCIRYEALGTPGFHYVDVVPPQATTFTLTGLQPSTRYRVWLLASNALG' \
              'DSGLADKGTQLPITTPGLHQPSGEPEDQLPTEPPSGPSGLPLLPVLFALGGLLLLSNASC' \
              'VGGVLWQRRLRRLAEGISEKTEAGSEEDRVRNEYEESQWTGERDTQSSTVSTTEAEPYYR' \
              'SLRDFSPQLPPTQEEVSYSRGFTGEDEDMAFPGHLYDEVERTYPPSGAWGPLYDEVQMGP' \
              'WDLHWPEDTYQDPRGIYDQVAGDLDTLEPDSLPFELRGHLV'
    # nicd = nephrin[1077:1242]
    # norm = NormScore([nicd])
    # print(norm.lc_norm_score())
    # print(nicd)
    # ns = tools_lc.display_lc(nicd, k, lca, lce)
    # print(ns)
    # motifs = tools_lc.count_lc_motifs(nicd, k, lca, lce)
    # print(motifs)
    # gfp = 'MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL' \
    #       'VTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLV' \
    #       'NRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLAD' \
    #       'HYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK'
    # #ns = tools_lc.display_lc(nicd+gfp, k, lca, lce)
    # #norm = NormScore([nicd+gfp])
    # #print(norm.lc_norm_score())
    # #print(ns)
    # fus = 'MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQ' \
    #       'SQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSS' \
    #       'SYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQY'
    # ns = tools_lc.display_lc(fus, k, lca, lce)
    # print(ns)
    # print(fus)
    # norm = NormScore([fus])
    # print(norm.lc_norm_score())




if __name__ == '__main__':
    main()