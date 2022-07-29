def drill(d):1
    if len(d) == 2 and isinstance(d[0], float):
        ans.append((d[0], d[1]))
        return
    for x in d:
         drill(x)

output = {}

for state in states["features"]:
   name = state["properties"]["NAME"]
   print(name)
   ans = []
   drill(state["geometry"]["coordinates"])
   output[name] = ans
   
