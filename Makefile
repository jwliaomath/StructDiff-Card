.PHONY: test demo site app

test:
	pytest

demo:
	structdiff compare examples/input/2HHB.pdb examples/input/1HHO.pdb -o examples/result

site:
	cd site && pnpm install --frozen-lockfile && pnpm build

app:
	streamlit run app.py