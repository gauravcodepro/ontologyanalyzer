def IPRNetworkPlot(tree_file, arg_IPR = None, arg_domain = None, arg_type = None):
    """
    _summary_
    A faster function to normalize the parent
    child relationship of the interpro domains and
    plot the undirected graph for the desired interpro
    domain. 

    Arguments:
        tree_file -- _description_
        a intepro domain parent child tree file
        for viewing the undirected acyclic network
        graph. The current release of the tree file
        can be found at: https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/ParentChildTreeFile.txt
    """
    if arg_type == "interpro" and arg_IPR:
        read_file=list(filter(None,[i.strip() for i in open(tree_file).readlines()]))
        domaindict={};
        for line in read_file:
            if line.startswith('IPR'):
                    doname = line.strip()
                    if line not in domaindict:
                        domaindict[line] = ""
                    continue
            domaindict[doname] += line.strip()
        normalize_domain = {}
        for k,v in domaindict.items():
            normalize_domain[k.split("::")[0]] = [i for i in ([i.replace("--", "") \
                                            if i.startswith("--") else i for i in \
                                                v.split("::")]) if i.startswith("IPR")]
        normalize_categories = {}
        for k,v in domaindict.items():
            normalize_categories[k.split("::")[1]] = list(filter(None,[i for i in \
                                                    ([i.replace("--", "") if i.startswith("--") \
                                                           else i for i in v.split("::")]) if \
                                                                         not i.startswith("IPR")]))
        return [v for k,v in normalize_domain.items() if k == arg_IPR]

    if arg_type == "domain" and arg_domain:
        read_file=list(filter(None,[i.strip() for i in open(tree_file).readlines()]))
        domaindict={};
        for line in read_file:
            if line.startswith('IPR'):
                    doname = line.strip()
                    if line not in domaindict:
                        domaindict[line] = ""
                    continue
            domaindict[doname] += line.strip()
        normalize_domain = {}
        for k,v in domaindict.items():
            normalize_domain[k.split("::")[0]] = [i for i in ([i.replace("--", "") \
                                            if i.startswith("--") else i for i in \
                                                v.split("::")]) if i.startswith("IPR")]
        normalize_categories = {}
        for k,v in domaindict.items():
            normalize_categories[k.split("::")[1]] = list(filter(None,[i for i in \
                                                    ([i.replace("--", "") if i.startswith("--") \
                                                           else i for i in v.split("::")]) if \
                                                                         not i.startswith("IPR")]))
        return [v for k,v in normalize_categories.items() if k == arg_domain]

def bacterialProteomeDomainAnalyzer(file, arg_type = None):
    import json 
    """sumary_line
    a rapid implementation of the bacterial domain
    analyzer, following the predictions of the bacterial
    domains from the interproscan. I implemented a mapped 
    dataframe approach to make it faster and iterable. it will 
    parse a nested to nested json from interpro for direct analysis
    and fecthing all the protein domains and the corresponding start
    and stop coordinates. Also you can make a direct ingestion to the
    database and it also provides a dataframe.
    Keyword arguments:
    argument -- file prediction by the interproscan
    Return: a systematic prediction of the domains in the sequences
    """
    if arg_type == "sequence":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        sequence = ''.join([i["sequence"] for i  in data["results"]])
        return sequence
    if arg_type == "interpro_normalize":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        interpro_normalize = pd.concat(list(map(lambda n: pd.DataFrame(n), \
                                [i["matches"] for i  in data["results"]])))
        return interpro_normalize
    if arg_type == "signature":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        signature = pd.concat(list(map(lambda n: pd.DataFrame(n),\
                                [i["matches"] for i  in data["results"]])))['signature']
        return signature
    if arg_type == "location":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        location = pd.concat(list(map(lambda n: pd.DataFrame(n),\
                                [i["matches"] for i  in data["results"]])))['locations']
        return location
    if arg_type == "evalue":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        evalue = pd.concat(list(map(lambda n: pd.DataFrame(n),\
                                [i["matches"] for i  in data["results"]])))['evalue']
        return evalue
    if arg_type == "score":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        score = pd.concat(list(map(lambda n: pd.DataFrame(n),\
                                [i["matches"] for i  in data["results"]])))['score']
        return score
    if arg_type == "modelac":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        modelac = pd.concat(list(map(lambda n: pd.DataFrame(n),\
                                [i["matches"] for i  in data["results"]])))['model-ac']
        return modelac
    if arg_type == "scope":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        scope = pd.concat(list(map(lambda n: pd.DataFrame(n),\
                               [i["matches"] for i  in data["results"]])))['scope']
        return scope
    if arg_type == "accession":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        accession = pd.concat(list(map(lambda n: pd.DataFrame(n),\
                               [i["matches"] for i  in data["results"]])))['accession']
        return accession
    if arg_type == "name":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        name = pd.concat(list(map(lambda n: pd.DataFrame(n),\
                              [i["matches"] for i  in data["results"]])))['name']
        return name
    if arg_type == "proteinClass":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        proteinClass = pd.concat(list(map(lambda n: pd.DataFrame(n),\
                              [i["matches"] for i  in data["results"]])))['ProteinClass']
        return proteinClass
    if arg_type == "graftPoint":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        graftPoint = pd.concat(list(map(lambda n: pd.DataFrame(n),\
                              [i["matches"] for i  in data["results"]])))['graftPoint']
        return graftPoint
    if arg_type == "goxRefs":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        goxRefs = pd.concat(list(map(lambda n: pd.DataFrame(n),\
                              [i["matches"] for i  in data["results"]])))['goxRefs']
        return goxRefs
    if arg_type == "prediction_locations":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        prediction_locations = pd.DataFrame.from_dict(pd.concat(list(map(lambda n:\
                            pd.DataFrame(n),[i["matches"] for i  in data["results"]])))\
                                                        ["signature"].apply(lambda n: n.values()))
        return prediction_locations
    if arg_type == "getdomains":
        interpro_file = file
        with open(interpro_file) as read:
            data = json.load(read)
        get_domains = list(list(filter(lambda n: n!=None and n!="",i)) for i in \
                        (list(list(filter(lambda n: not isinstance(n,dict),i)) \
                                for i in ([list(i) for i in pd.DataFrame.from_dict \
                                            (pd.concat(list(map(lambda n: pd.DataFrame \
                                                (n),[i["matches"] for i  in data["results"]]))) \
                              ["signature"].apply(lambda n: n.values()))["signature"].to_list()]))))
        return get_domains
