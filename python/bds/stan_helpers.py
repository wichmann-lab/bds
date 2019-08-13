import pystan
import pickle
from hashlib import md5

def StanModel_cache(model_code, model_name=None, **kwargs):
    """Use just as you would `stan`

    taken from:
    https://pystan.readthedocs.io/en/latest/avoiding_recompilation.html
    """
    code_hash = md5(model_code.encode('ascii')).hexdigest()
    if model_name is None:
        cache_fn = 'cached-model-{}.pkl'.format(code_hash)
    else:
        cache_fn = 'cached-{}-{}.pkl'.format(model_name, code_hash)
    try:
        sm = pickle.load(open(cache_fn, 'rb'))
    except:
        sm = pystan.StanModel(model_name = model_name, model_code=model_code)
        with open(cache_fn, 'wb') as f:
            pickle.dump(sm, f)
    else:
        print("Using cached StanModel")
    return sm

def fit_stan_model(model_dir, model_name, data, init='random'):
  with open(model_dir + model_name + '.stan', 'r') as model_file:
    model_code = model_file.read()
  model = StanModel_cache(model_code, model_name=model_name)

  fit = model.sampling(data=data,
                       iter=2000,
                       chains=4,
                       control={'adapt_delta': 0.99},
                       init=init)

  return fit
