import metabo_engineer as me

m, _ = me.cargar_modelo_universal('iML1515')

discrepancias = ['b0040', 'b0312', 'b0655', 'b0831', 'b0864', 'b1022', 
                 'b2361', 'b2788', 'b3781', 'b4006', 'b4197']

print(f'Genes con discrepancia: {len(discrepancias)}')
print()
print(f"{'Gen':10s} | {'Reacciones':10s} | {'Estado':20s}")
print('-' * 50)

for g_id in discrepancias:
    if g_id in m.genes:
        gen = m.genes.get_by_id(g_id)
        num_rxn = len(gen.reactions)
        tiene = 'SÍ' if num_rxn > 0 else 'NO (gen huérfano)'
        print(f'{g_id:10s} | {num_rxn:10d} | {tiene:20s}')
    else:
        print(f'{g_id:10s} | NO EXISTE en el modelo')