<template>
  <div class="top">
    <table>
      <tr>
        <td>
          <img :src="'static/'+top.logo" />
        </td>
        <td>
          <span class="top-title">NSP</span><br><span class="top-intro">Nucleic acids Structure Prediction</span>
        </td>
      </tr>
    </table>
  </div>
</template>

<script>
  export default {
    data() {
      return {
        top: []
      }
    },
    created() {
      this.fetch_data()
    },
    // watch: {
    //   '$route': 'fetch_data'
    // },
    methods: {
      fetch_data() {
        var v = this
        this.$http.get('static/top.json').then((response) => {
          v.top = response.body
        }, (response) => {
        })
      }
    },
  }
</script>

<style>
.top {
  text-align: left;
  height: 50px;
  padding: 10px 0px;
}
.top img {
  height: 50px;
}
.top-title {
    font-size: 110%;
}
.top-intro {
    font-size: 80%;
}
</style>